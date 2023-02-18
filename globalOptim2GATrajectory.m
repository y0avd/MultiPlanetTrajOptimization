clear, clc, close all;

GAoptions = optimoptions('ga','MaxStallGenerations',10,...
    'PopulationSize',50);
GAoptions.EliteCount = 20;
GAoptions.Display = 'iter';
GAoptions.FunctionTolerance = 1e-12;
GAoptions.ConstraintTolerance = 1e-6;
%GAoptions.CreationFcn = 'gacreationlinearfeasible';
GAoptions.FitnessScalingFcn = 'fitscalingrank';
GAoptions.UseParallel = true;

dt = datetime('01-01-2031','InputFormat','dd-MM-yyyy');
launchDateRange = [juliandate(dt)-yr2day(1),juliandate(dt)+yr2day(1)];
sequence = [3,2,3,7];
uranus_vinf = 18;

lb = [0, 0.2, 0.2, 0.2];
ub = [1, 1, 1, 1];
nvars = 4;

maxTOF = yr2day(15);
minTOF = yr2day(10);

A = [0,getHohmannTOF(sequence(1),sequence(2)),...
    getHohmannTOF(sequence(2),sequence(3)),...
    getHohmannTOF(sequence(3),sequence(4));0,...
    -getHohmannTOF(sequence(1),sequence(2)),...
    -getHohmannTOF(sequence(2),sequence(3)),...
    -getHohmannTOF(sequence(3),sequence(4))];

b = [maxTOF,-minTOF];

exitflag = -1;

while exitflag < 0
    [x,fval,exitflag,output,population,scores] =...
    ga(@(input)twoGA_Trajectory(input,launchDateRange,sequence,2),...
    nvars,A,b,[],[],lb,ub,...
    @(input)nonLCons(input,launchDateRange,sequence,uranus_vinf),GAoptions);

    [c,ceq] = nonLCons(x,launchDateRange, sequence, uranus_vinf);
end

printOrbitInfo(x,launchDateRange,sequence)
plotTraj(x,launchDateRange,sequence)