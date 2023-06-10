clc
clear all
close all


%% inputs

% network
nodes = 10;
load('network.mat'); %update the path according to your PC

% parametrs of optimization problem (Target tracking via least-squares)
trgt = 3;
dim = 2*trgt; %dimension of the optimization variable

load('Problem.mat'); %%update the path according to your PC

% Algorithm
ratio = [0.1:0.1:2];  % algorithm parameter is defined as ratio*(1/L) 

%% Matrix formulation

L1 = L_mat;
M1 = mu_mat;

M = [-2*(L1*M1), L1 + M1 ;      %for strongly convex function
      L1 + M1, -2*eye(nodes)];

L =  (1/nodes)*sum(sum(L_mat));
mu = (1/nodes)*sum(sum(mu_mat));

%%% finding P with min rho for DEXTRA algorithm %%%


for i=1:length(ratio)
    i
    
    %DEXTRA
    DExeta = ratio(i)/L;
    M_1=zeros(nodes,nodes);
    M_2=eye(nodes);
    M_3=-0.5*(eye(nodes)+W);
    M_4= (eye(nodes)+W);
    M_5= DExeta*eye(nodes);
    F=[-1*ones(1,nodes) 1*ones(1,nodes) DExeta*ones(1,nodes)];
    G=zeros(1,nodes);
   [DExc_mat(i).P, DExcrho(i),DExclmbda1(i), DExc_cvx(i).s, DExc_t(i)] = fsblty_of_algo(nodes,dim,M,M_1,M_2,M_3,M_4,M_5,F,G); 
                                                                 % please first update the function code for strongly convex case
                                                                 % by uncommenting the corresponding part of the function
                                                             
    %DGT
    DGTeta = ratio(i)/L;
    M_1 = W;
    M_2 = W;
    M_3 = zeros(nodes);
    M_4 = W;
    M_5 = DGTeta*eye(nodes);
    F   = [zeros(1,nodes) 1*ones(1,nodes) DGTeta*ones(1,nodes)];
    G   = zeros(1,nodes);
    [DGTc_mat(i).P,DGTcrho(i),DGTclmbda1(i), DGTc_cvx(i).s, DGTc_t(i)] = fsblty_of_algo(nodes,dim,M,M_1,M_2,M_3,M_4,M_5,F,G);
 
end


%saving factor

DExrho = DExcrho;
DExfactrfnl = ((DExrho.*DExrho)./((1-DExrho).^2)).*DExc_t;
DExtuned_eta = ratio(find(DExfactrfnl==min(DExfactrfnl)))/L;


DGTrho = DGTcrho;
DGTfactrfnl = ((DGTrho.*DGTrho)./((1-DGTrho).^2)).*DGTc_t;
DGTtuned_eta = ratio(find(DGTfactrfnl==min(DGTfactrfnl)))/L;

fname = sprintf('Tuned_hyp_par');
save(fname) 


% for making rho less than 1



% %% plotting figure
% 
% str = '#D95319';
% color1 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
% 
% 
% %plotting rho for all algorithms
% figure(1)
% plot(ratio,DExrho,"Color",'r','linewidth',2);
% xlabel('stepsize(\eta)*L','FontSize',12,'FontWeight','bold')
% ylabel('\rho','FontSize',12,'FontWeight','bold')
% hold on
% plot(ratio,DNsrho,"Color",'g','linewidth',2);
% plot(ratio,DPGPrho,"Color",'c','linewidth',2);
% plot(ratio,DGErho,"Color",'m','linewidth',2);
% plot(ratio,DATrho,"Color",'y','linewidth',2);
% plot(ratio,DGTrho,"Color",'b','linewidth',2);
% plot(ratio,DooGTrho,"Color",color1,'linewidth',2);
% plot(ratio,oGGTrho,"Color",'k','linewidth',2);
% legend('DEXTRA','DNIDS','DPGPDA','DGE','DATCGT','DGT','DOOGT','oGGT','Fontsize',10,'Fontweight','bold');
% grid on




