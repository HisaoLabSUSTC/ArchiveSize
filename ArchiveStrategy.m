classdef ArchiveStrategy < handle
    properties
        % 1 - standard strategy
        % 2 - lazy periodical strategy
        % 3 - last X-generation strategy
        strategy

        K % size of archive
        T % update interval
        cnt % number of solutions in the archive
        preCnt % number of examined solutions in the archive
        name % name of the strategy
        dataName % name of the sequence
        runid % run of the sequence
        saveName % name of the save file 
        % To save the storage, the index of solution is recored. postive 
        % for the offspring, and negative for the population
        arcIndex = [];
        archive = [];
        % Removal methods have two options: "List" means using a linear  
        % list and "ENS" means using T-ENS method.
        removalMethod = "List"; 
        % Truncation methods have two options: "DSS" means using greedy 
        % distance-based inclusion and "CDR" means using greedy crowding
        % distance-based removal.
        truncationMethod = 'DSS';
        folder = "";
    end
    
    methods
        function obj = ArchiveStrategy(strategy, K, T, folder)
            obj.strategy = strategy;
            obj.K = K;
            obj.T = T;
            obj.folder = folder;
            obj.name = sprintf("S%d", obj.strategy);
            if strategy==1 || strategy==2
                obj.name = obj.name + sprintf("_K%d_T%d", K, T); 
            elseif obj.strategy==3
                obj.name = obj.name + sprintf("_K%d", K);
            end
        end
        function update(obj, seq)
            start = tic; % timer (limied in 1 hour)
            gen = size(seq,2);
            [N, M] = size(seq(1).population);
            % num1 - the number of solutions in the archive before dominted 
            %        solution removal; 
            % num2 - the number solutions in the archive after truncation;
            num1 = zeros(1,gen);
            num2 = zeros(1,gen);
            removalTime = zeros(1,gen);
            truncationTime = zeros(1,gen);

            % standard strategy and lazy periodical strategy
            if obj.strategy == 1 || obj.strategy == 2
                % pre-allocate memory
                memory = obj.K+obj.T*N;
                obj.archive = zeros(memory, M);
                obj.cnt=0;
                obj.preCnt=0;
                % A_1=P_1
                obj.archive(1:N,:)=seq(1).population;
                obj.cnt=N;
                obj.arcIndex = -[1:N];
                num1(1)=N;  
                num2(1)=N; 
                removalTime(1)=0;
                truncationTime(1)=0;

                solIndex = 0; 
                for g=2:gen
                    % A_g=A_g-1+O_g-1
                    O = seq(g-1).offspring;
                    n = size(O, 1); % n maybe smaller than N (e.g., n=90 and N = 91)
                    obj.archive(obj.cnt+1: obj.cnt+n,:) = O;
                    obj.cnt = obj.cnt + n;
                    obj.arcIndex = [obj.arcIndex, solIndex+1:solIndex+n];
                    solIndex = solIndex + n;
                    num1(g)=obj.cnt;
                    % standard strategy
                    if obj.strategy == 1
                        % remove dominated solutions
                        removalTime(g)=obj.removeD();
                        % truncate
                        if obj.cnt>obj.K
                            truncationTime(g)=obj.truncate();
                        end
                    % lazy periodical strategy
                    elseif obj.strategy == 2
                        if mod((gen-g),obj.T)==0
                            % remove dominated solutions
                            if obj.cnt>obj.K || g==gen
                                removalTime(g)=obj.removeD();
                            end
                            % truncate
                            if obj.cnt>obj.K
                                truncationTime(g)=obj.truncate(); 
                            end
                        end
                    end
                    num2(g) = obj.cnt;
                    if sum(truncationTime)+sum(removalTime) > 3600
                        obj.saveEmpty();
                        return;
                    end
                end
            end
            % last X-generation strategy
            if obj.strategy == 3 
                X = min(floor(obj.K/N),gen);
                memory = X * N;
                % pre-allocate memory
                obj.archive = zeros(memory,M);
                obj.cnt=0;
                obj.preCnt=0;
                solIndex_O=0;
                solIndex_P=0;
                for g=1:gen-X
                    solIndex_O = solIndex_O + size(seq(g).offspring,1);
                    solIndex_P = solIndex_P + size(seq(g).population,1);
                end
               
                for g=(gen-X+1):gen
                    if g==gen-X+1 % initialize the archive
                        P=seq(g).population;
                        n = size(P,1);
                        obj.archive(1:n,:) = P;
                        obj.cnt = n;
                        obj.arcIndex=-[solIndex_P+1:solIndex_P+n];
                    else
                        O = seq(g-1).offspring;
                        n = size(O, 1);
                        obj.archive(obj.cnt+1: obj.cnt+n, :)=O;
                        obj.cnt = obj.cnt + n;
                        obj.arcIndex=[obj.arcIndex, solIndex_O+1:solIndex_O+n];
                        solIndex_O = solIndex_O + n;
                    end
                    num1(g)=obj.cnt;
                    if g==gen
                        removalTime(g)=obj.removeD();   
                    end
                    num2(g)=obj.cnt;
                    % note
                    if sum(truncationTime)+sum(removalTime) > 3600
                        obj.saveEmpty();
                        return;
                    end
                end
            end
            % save the archive for final selection
            arcIndex = obj.arcIndex;
            runtime = sum(truncationTime)+sum(removalTime);
            save(obj.saveName, ...
                'arcIndex', 'memory', 'num1', 'num2', 'truncationTime','removalTime','runtime');
        end
        function res = normalize(~, val, PF) % nomrmalize the obecjtive values by PF
            N = size(val,1);
            fmin   = min(PF, [], 1);
            fmax   = max(PF,[],1);
            res = (val-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);
        end

        function time = removeD(obj) 

             if obj.removalMethod == "ENS"
                 arc = obj.archive(1:obj.cnt, :);
                 st=tic;
                 [FrontNo,~] = NDSort(arc, 1);   
                 arc = arc(FrontNo==1,:);
                 time = toc(st);
                 obj.cnt = size(arc,1);
                 obj.preCnt = obj.cnt;
                 obj.archive(1:obj.cnt,:) = arc;
                 obj.arcIndex = obj.arcIndex(FrontNo==1); 
             elseif obj.removalMethod == "List"
                % only nondominated solutions in S are added to A
                A = obj.archive(1:obj.preCnt, :);
                Aindex = obj.arcIndex(1:obj.preCnt);
                S = obj.archive(obj.preCnt+1:obj.cnt,:);
                Sindex = obj.arcIndex(obj.preCnt+1:obj.cnt);

                st=tic;
                [FrontNo, ~]  = NDSort(S,1);
                S = S(FrontNo==1,:);
                Sindex = Sindex(FrontNo==1);
                p = size(S,1);
                k = obj.preCnt;

                if k>0 && p>0
                    IA = true(1,k);
                    IS = true(1,p);
                    for i=1:p
                        value = A-repmat(S(i,:), k, 1);
                        % remove dominated solutions in the archive
                        IA(all(value>=0, 2) & (~all(value==0, 2)))=false;
                        % remove the dominated solution in S
                        if any(all(value<=0, 2),1)
                            IS(i)=false;
                        end
                    end
                    S = S(IS,:);
                    Sindex = Sindex(IS);
                    A = A(IA, :);
                    Aindex = Aindex(IA);
                end
                time = toc(st);
                
                p = size(S,1);
                k = size(A,1);
                obj.cnt = k + p;
                obj.preCnt = obj.cnt;
                obj.archive(1:k+p,:)=[A;S];
                obj.arcIndex=[Aindex,Sindex];
             else
                 error('removal method is not implemented.')
             end
        end

        function time = truncate(obj)
             arc= obj.archive(1:obj.cnt,:);
             st = tic;
             if strcmp(obj.truncationMethod, 'DSS')
                 [~,index] = DSS(obj.normalize(arc, arc), obj.K);
             elseif strcmp(obj.truncationMethod, 'CDR')
                 index = CDR(obj.normalize(arc, arc), obj.K);
             else
                 error('truncation method is not implemented.')
             end
             arc=arc(index,:);
             time = toc(st);
             obj.cnt=size(arc,1);
             obj.preCnt = obj.cnt;
             obj.archive(1:obj.cnt,:) = arc;
             obj.arcIndex = obj.arcIndex(index); 
        end

        function saveEmpty(obj)
            runtime = -1;
            save(obj.saveName, "runtime");
        end

        function done = setName(obj, dataName, runid)
            obj.dataName = dataName;
            obj.runid = runid;
            if (~exist(sprintf('%s/%s/',obj.folder, obj.dataName),'dir'))
                mkdir(sprintf('%s/%s/',obj.folder, obj.dataName))
            end
            done =  (~exist(sprintf('%s/%s/%s_%d.mat',obj.folder,dataName,obj.name, runid),'file')==0);
            obj.saveName = sprintf('%s/%s/%s_%d.mat',obj.folder,  obj.dataName, obj.name, obj.runid);
        end
    end
end

