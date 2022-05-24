classdef NSGAIII_ARC < ALGORITHM
% <multi/many> <real/binary/permutation> <constrained/none>
% Nondominated sorting genetic algorithm III

%------------------------------- Reference --------------------------------
% K. Deb and H. Jain, An evolutionary many-objective optimization algorithm
% using reference-point based non-dominated sorting approach, part I:
% Solving problems with box constraints, IEEE Transactions on Evolutionary
% Computation, 2014, 18(4): 577-601.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate the reference points and random population
            startTime = tic;
            [Z, Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population    = Problem.Initialization();
            Zmin          = min(Population(all(Population.cons<=0,2)).objs,[],1);
            index = Algorithm.ParameterSet(1);
            %% Optimization
            %% Optimization
            n = 1;
            rec(n).time = toc(startTime);
            rec(n).population = Population.objs;
            rec(n).FE = Problem.FE;
            while Algorithm.NotTerminated(Population)
                startTime = tic;
                MatingPool = TournamentSelection(2,Problem.N,sum(max(0,Population.cons),2));
                Offspring  = OperatorGA(Population(MatingPool));
                Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
                Population = EnvironmentalSelection([Population,Offspring],Problem.N,Z,Zmin);
                n = n + 1;
                rec(n).time = toc(startTime);
                rec(n).population = Population.objs;
                rec(n-1).offspring = Offspring.objs;
                rec(n).FE = Problem.FE;
            end
            save(sprintf('/home/vmuser/Dataset/ArchiveSize/Data/NSGAIII_%s_N%d_M%d_%d.mat', class(Problem), Problem.N, Problem.M, index),"rec");
        end
    end
end