algos = ["NSGAII", "NSGAIII", "MOEAD_PBI"];
problems = ["DTLZ1","DTLZ2","DTLZ3","DTLZ4","MinusDTLZ1", "MinusDTLZ2", "MinusDTLZ3", "MinusDTLZ4"];
objs = [3,5,8];%[3,5,8];
n=0;
PF =  containers.Map;
Gens = [400,600,750; 250,350,500; 1000,1000,1000; 600,1000,1250;400,600,750; 250,350,500; 1000,1000,1000; 600,1000,1250];
arcs = [1,1;2,1]; %1,1; 2,1;
method = ["DSS", "LazyHV"];
folder = "/home/vmuser/Dataset/ArchiveSize/";
%folder="/home/vmuser/Nutstore Files/2022/project/TEVC_ArchiveSize/PlatEMO-master/PlatEMO/tmp";
savefolder = sprintf("%s/FinalArchive/",folder);
%Get the PF for each experiment
for b=1:length(problems)
    for c=1:length(objs)
       % sprintf('PF(\"%s_%d\")=%s(''M'',%d);',problems(b), objs(c), problems(b),objs(c))
        eval(sprintf('PF(\"%s_%d\")=%s(''M'',%d).optimum;',problems(b), objs(c), problems(b),objs(c)));
    end
end
clear task
for b=1:length(problems)
    for c=1:length(objs)
        for e=1:size(arcs,1)
            if ~(ismember(num2str(objs(c))+problems(b),...
                ["8DTLZ3"]))
                continue;
            end

            if objs(c)==3
                pop = 91;
            elseif objs(c)==5
                pop = 210;
            else
                pop = 156;
            end
            for arcsiz = [1,2,5,10,20,50,100,200,500,1000,2000]*pop
                clear task
                n=0;
                for a=1:length(algos)
                    for d=1:21   
                        name = sprintf("%s_%s_M%d",algos(a), problems(b),objs(c));
                        gens=Gens(b,c);
                        S = arcs(e,1); T=arcs(e,2);
                        arcst = ArchiveStrategy(S,arcsiz,T, savefolder);
                        arcst.setName(name, d);
                        name2 = sprintf("%s_%d",arcst.name,d);
                        if (~exist(sprintf('%s/Nadir/%s',folder, name),'dir'))
                            mkdir(sprintf('%s/Nadir/%s/',folder, name))
                        end
                        n=n+1;
                        task(n).algo = algos(a);
                        task(n).problem = problems(b);
                        task(n).obj = objs(c);
                        task(n).id = d;
                        task(n).pop = pop;
                        task(n).loadPath = arcst.saveName;
                        task(n).savePath=sprintf('%s/Nadir/%s/%s.mat',folder,name,name2);
                    end
                end
                parfor i = 1:63
                    index = i;
                    algo = task(index).algo;
                    problem = task(index).problem;
                    obj = task(index).obj;
                    pop = task(index).pop;
                    id = task(index).id;
                    loadPath = task(index).loadPath;
                    savePath = task(index).savePath;
                    name = sprintf("%s_%s_N%d_M%d_%d",algo, problem, pop, obj, id);
                    seq = load(sprintf("%s/Data/%s","/home/vmuser/Dataset/ArchiveSize",name)).rec;
                    result = load(loadPath);
                    selection(seq, result,PF(sprintf("%s_%d",problem,obj)), savePath);
                    sprintf("Work%d is done.", i);
                end
            end
        end   
    end
end

function selection(seq, result, PF, path)
    population = [];
    offspring = [];
    N = size(seq(1).population,1);
    for i=1:size(seq,2)
        population = [population;seq(i).population];
        offspring = [offspring;seq(i).offspring];
    end
    if result.runtime==-1 % exceed the time limit
        runtime=-1;
        save(path,"runtime"); 
        return;
    end

    arcIndex = result.arcIndex;
    truncationTime = sum(result.truncationTime);
    removalTime = sum(result.removalTime);
    arc = [offspring(arcIndex(arcIndex>0),:); population(-arcIndex(arcIndex<0),:)];

    tic;
    nadir1 = max(arc,[],1);
    solnum1 = size(arc,1);
    [~,index] = LazyHVSelection(normalize(arc, arc), N, ones(1,size(arc,2))*1.2);
    selectionTime = toc;
    runtime = removalTime+truncationTime+selectionTime;
    if runtime>3600
        runtime=-1;
        save(path,"runtime");
        return;
    end
    res = arc(index,:);
    normalizedRes = normalize(res, PF);
    hv1 = HV(normalizedRes,1.2);

    meanVal = repmat(mean(arc,1), size(arc,1),1);
    stdVal = repmat(std(arc,1), size(arc,1),1);
    arc(any((arc-meanVal)./stdVal>6,2),:)=[];
    nadir2 = max(arc,[],1);
    solnum2 = size(arc,1);
    
    [~,index] = LazyHVSelection(normalize(arc, arc), N, ones(1,size(arc,2))*1.2);
    res = arc(index,:);
    normalizedRes = normalize(res, PF);
    hv2 = HV(normalizedRes,1.2);   
    save(path,"runtime","nadir1","res","nadir2", "solnum1", "solnum2", "hv1","hv2");
end
function res = normalize(val, PF) % nomrmalize the obecjtive values by PF
    N = size(val,1);
    fmin   = min(PF, [], 1);
    fmax   = max(PF,[],1);
    res = (val-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);
end
