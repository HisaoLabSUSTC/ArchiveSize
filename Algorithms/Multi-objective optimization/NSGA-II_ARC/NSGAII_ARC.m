classdef NSGAII_ARC < ALGORITHM
% <multi> <real/binary/permutation> <constrained/none>
% Nondominated sorting genetic algorithm II

%------------------------------- Reference --------------------------------
% K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, A fast and elitist
% multiobjective genetic algorithm: NSGA-II, IEEE Transactions on
% Evolutionary Computation, 2002, 6(2): 182-197.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm, Problem)
            %% Generate random population
            startTime = tic;
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);
            index = Algorithm.ParameterSet(1);
            %% Optimization
            n = 1;
            rec(n).time = toc(startTime);
            rec(n).population = Population.objs;
            rec(n).FE = Problem.FE;
            while Algorithm.NotTerminated(Population)
                startTime = tic;
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorGA(Population(MatingPool));
                [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Problem.N);
                n = n + 1;
                rec(n).time = toc(startTime);
                rec(n).population = Population.objs;
                rec(n-1).offspring = Offspring.objs;
                rec(n).FE = Problem.FE;
            end
            save(sprintf('/home/vmuser/Dataset/ArchiveSize/Data/NSGAII_%s_N%d_M%d_%d.mat', class(Problem), Problem.N, Problem.M, index),"rec");
            %save(sprintf('ArcExp/PaperResult/Example/NSGAII_%s_N%d_M%d_%d.mat', class(Problem), Problem.N, Problem.M, index),"rec");
        end
    end
end