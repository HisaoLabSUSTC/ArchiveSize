paras = [3,91; 5,210; 8,156]; % objective
problems = {@DTLZ1, @DTLZ2, @DTLZ3, @DTLZ4, @MinusDTLZ1, @MinusDTLZ2, @MinusDTLZ3, @MinusDTLZ4};
Gens = [400,600,750; 250,350,500; 1000,1000,1000; 600,1000,1250;400,600,750; 250,350,500; 1000,1000,1000; 600,1000,1250];
algos = {@NSGAII_ARC, @NSGAIII_ARC, @MOEAD_ARC};

parfor i = 1:21
    for a = 1:length(algos)
        for b = 1:length(problems)
            for c = 1:size(paras,1)
                M = paras(c,1);
                N = paras(c,2);
                if a<=2 && mod(N,2)==1
                    % the offspring of NSGAII and NSGAIII is multiple of 2
                    maxFE = (N-1)*Gens(b,c)+1; 
                else
                    maxFE = paras(c,2)*Gens(b,c);
                end
                platemo('algorithm',{algos{a},i},'problem',problems{b}, 'N', N,'M',M, 'maxFE', maxFE);
            end
        end
    end
end