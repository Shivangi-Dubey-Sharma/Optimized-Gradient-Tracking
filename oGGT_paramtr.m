clc
clear all
close all


nodes = 10; 
load('network.mat'); %update the path according to your PC
% parametrs of optimization problem (Target tracking via least-squares)
trgt = 3;
dim = 2*trgt; %dimension of the optimization variable

load('Problem.mat'); %%update the path according to your PC
load('Tuned_hyp_par.mat'); %%update the path according to your PC

%% Matrix formulation

L1 = L_mat;
M1 = mu_mat;

M = [-2*(L1*M1), L1 + M1 ;      %for strongly convex function
      L1 + M1, -2*eye(nodes)];

L =  (1/nodes)*sum(sum(L_mat));
mu = (1/nodes)*sum(sum(mu_mat));


fac=1; %to vary the range of parameters
alpha = optimizableVariable('al',[0,1]*fac,'Type','real');
beta  = optimizableVariable('be',[0,1]*fac,'Type','real');
gama  = optimizableVariable('ga',[0,1]*fac,'Type','real');
delta = optimizableVariable('de',[0,1]*fac,'Type','real');


M_1= W;
M_2= @(x)x.al*eye(nodes)+x.be*W;
M_3= zeros(nodes);
M_4= W;
M_5= @(x)x.ga*eye(nodes)+x.de*W;
F=@(x)[zeros(1,nodes) 1*ones(1,nodes) ones(1,nodes)*M_5(x)];
G=zeros(1,nodes);

fun = @(x)(fsblty_of_algo_bysopt(nodes,dim,M,M_1,M_2(x),M_3,M_4,M_5(x),F(x),G,[0.6 1]).rho);
results = bayesopt(fun,[alpha,beta,gama,delta],'MaxObjectiveEvaluations',1000)
