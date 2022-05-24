algos = ["NSGAII", "NSGAIII", "MOEAD_PBI"];
problems =["DTLZ1","DTLZ2","DTLZ3","DTLZ4","MinusDTLZ1", "MinusDTLZ2", "MinusDTLZ3", "MinusDTLZ4"];
objs = [3,5,8];%[3,5,8]; % all should be preserve
PF =  containers.Map;
Gens = [400,600,750; 250,350,500; 1000,1000,1000; 600,1000,1250;400,600,750; 250,350,500; 1000,1000,1000; 600,1000,1250];
arcs = [1,1; 2,1; 2,2; 2,5; 2,10;2,20; 3,1];%[1,1; 2,1; 2,2; 2,5; 2,10;2,20; 3,1];
folder = "/home/vmuser/Dataset/ArchiveSize/";
savefolder = sprintf("%s/FinalArchive/",folder);
%Get the PF for each experiment
for b=1:length(problems)
    for c=1:length(objs)
       % sprintf('PF(\"%s_%d\")=%s(''M'',%d);',problems(b), objs(c), problems(b),objs(c))
        eval(sprintf('PF(\"%s_%d\")=%s(''M'',%d).optimum;',problems(b), objs(c), problems(b),objs(c)));
    end
end



for b=1:length(problems)
    for c=1:length(objs)
        for e=1:size(arcs,1)
%             if (ismember(num2str(objs(c))+problems(b),...
%                 ["3MinusDTLZ1", "5DTLZ3", "8MinusDTLZ2"]))
%                 continue
%             end
            if ~(ismember(num2str(objs(c))+problems(b),...
                ["5MinusDTLZ3", "8MinusDTLZ3", "3MinusDTLZ4", "5MinusDTLZ4","8MinusDTLZ4"]))
                continue
            end
            if ismember(num2str(objs(c))+problems(b),...
                ["5MinusDTLZ3"]) && e<=4
                continue
            end
%             if ismember(num2str(objs(c))+problems(b),["5DTLZ3"]) && e==1
%                 continue
%             end
            if objs(c)==3
                pop = 91;
            elseif objs(c)==5
                pop = 210;
            else
                pop = 156;
            end
            for arcsiz = [1,2,5,10,20,50,100,200,500,1000,2000]*pop
%                 if ~(ismember(num2str(objs(c))+problems(b),...
%                  ["8MinusDTLZ2"]) && arcsiz==500*pop && e==1)
%                     continue
%                 end
                clear task
                n=0;
                for a=1:length(algos)
                    for d=1:21     
                        name = sprintf("%s_%s_M%d",algos(a), problems(b),objs(c));
                        gens=Gens(b,c);
                        S = arcs(e,1); T=arcs(e,2);
                        arcst = ArchiveStrategy(S, arcsiz, T, savefolder);
                        done = arcst.setName(name, d);
                        n=n+1;
                        task(n).algo = algos(a);
                        task(n).problem = problems(b);
                        task(n).obj = objs(c);
                        task(n).id = d;
                        task(n).pop = pop;
                        task(n).para = {S, arcsiz, T,  savefolder};
                    end
                end
                parfor i = 1:63
                    index = i;
                    algo = task(index).algo;
                    problem = task(index).problem;
                    obj = task(index).obj;
                    pop = task(index).pop;
                    id = task(index).id;
                    name = sprintf("%s_%s_N%d_M%d_%d",algo, problem, pop, obj, id);
                    data = load(sprintf("%s/Data/%s",folder,name));
                    para = task(index).para;
                    arcst = ArchiveStrategy(para{1}, para{2}, para{3}, para{4});
                    arcst.setName(sprintf("%s_%s_M%d",algo, problem, obj), id);
                    arcst.name
                    arcst.update(data.rec);
                    sprintf("Work%d is done.", i);
                end
            end
        end   
    end
end
