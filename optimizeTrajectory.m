clear, clc, close all;

dt = datetime('01-02-2031','InputFormat','dd-MM-yyyy');
launchDateRange = [juliandate(dt)-15,juliandate(dt)+15];
sequence = [3,4,5,7];
uranus_vinf = 6;

IC = [0.85    0.3683    0.6106    0.4579];
lb = IC - 0.1;
ub = IC + 0.1;

maxTOF = yr2day(15);
minTOF = yr2day(10);

A = [0,getHohmannTOF(sequence(1),sequence(2)),...
    getHohmannTOF(sequence(2),sequence(3)),...
    getHohmannTOF(sequence(3),sequence(4));0,...
    -getHohmannTOF(sequence(1),sequence(2)),...
    -getHohmannTOF(sequence(2),sequence(3)),...
    -getHohmannTOF(sequence(3),sequence(4))];

b = [maxTOF,-minTOF];

options = optimoptions("fmincon",'Display','iter','Algorithm','sqp');

x = fmincon(@(input)twoGA_Trajectory(input,launchDateRange,sequence,2),...
    IC,A,b,[],[],lb,ub,...
    @(input)nonLCons(input,launchDateRange,sequence,uranus_vinf),options);

printOrbitInfo(x,launchDateRange,sequence)
plotTraj(x,launchDateRange,sequence) 