%SIMPLE GENETIC ALGORITHM IMPLEMENTATION
clear all
close all

%set initial simulink parameters 
%simscapeModelParameters_16_08_tuning

pop = 10;
options = optimoptions('ga','Display','iter','PopulationSize',pop,'TolFun',1e-6,'useparallel',false,'UseVectorized',false);
% initPop = repmat([100 0.5 10], 3, 1);  %  individuals start near [1,1,1]
% options = optimoptions(options,'InitialPopulationMatrix',initPop);

A = [];
b = [];
Aeq = [];
beq = [];
lb = [1e-9 1e-9 1];
ub = [100 100 100];
nonlcon = [];


Tsim = 10
load_system('simplePatientModelTuning_01_11_23');
x = ga(@optimize_vent_control_param,3,A,b,Aeq,beq,lb,ub,nonlcon,options);
[x(1) x(2) x(3)]