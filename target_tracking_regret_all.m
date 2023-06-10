clc
clear all
close all

seed = 1.1;
rng(seed);

nodes = 10;
load('network.mat'); %update the path according to your PC

% parametrs of optimization problem (Target tracking via least-squares)
trgt = 3;
dim = 2*trgt; %dimension of the optimization variable

load('Problem.mat'); %%update the path according to your PC
load('Tuned_hyp_par.mat'); %%update the path according to your PC

max_t = 10000;

amp_min = [0.4 0.4 0.4]; 
amp_max = [0.8 0.8 0.8];

f = [100 100 100]; %input('Please enter the sampling frequency of the system in Hz : ');

prdcty = [50 50 50]; %input('Please enter the preodicity of each signal in sec in a row vector : ');


%% Utilizing stepsize which gives minimum rho

DExeta1 = DExtuned_eta; 
DNseta1 = DNstuned_eta;
DPGPc1 = DPGPtuned_eta; 
DGEeta1 = DGEtuned_eta; 
DATCGTeta1 = DATtuned_eta;  
DGTeta1 =  DGTtuned_eta; 
DooGTeta1 =  DooGTtuned_eta;

% for oGGT
alpha =  3.0525/L; 
beta =   20.1284/L; 
gamma =  0.7934/L; 
delta =  4.1387/L; 


%% data generation
x_data = zeros(2*trgt,max_t+1);
y_data = zeros(4*nodes,max_t+1);    

[x_star,y_d] = trgt_trckng_dataGen(nodes,trgt,amp_min,amp_max,prdcty,f,1:max_t+1,C_data);
x_data = x_star;
y_data = y_d;


%% runnig algorithm for each timestep
ts=0;

%x_zero = rand(2*trgt*nodes,1);  %can have random initialization also
x_zero = zeros(2*trgt*nodes,1);  % Just to replicate the results

nab = [];
for i=1:nodes
    nab = [nab; C_data((i-1)*6+1:i*6,:)'*(C_data((i-1)*6+1:i*6,:)*x_zero((i-1)*2*trgt+1:i*2*trgt,1) - y_data((i-1)*6+1:i*6,1))];
end
nab_zero = nab;


%Defining_initial_vectors
x_kp1 = x_zero - DExeta1*nab_zero;
x_k = x_zero;
DExz_kp1 = [x_k ; x_kp1; nab_zero];                %inital vector for DEXTRA

x_kp1 = x_zero - DNseta1*nab_zero;
x_k = x_zero;
DNsz_kp1 = [x_k ; x_kp1; nab_zero];                %inital vector for DNIDS


x_kp1 = x_zero - 0.5*(DPGPc1)*inv(kron(diag(d),eye(2*trgt)))*nab_zero;
x_k = x_zero;
DPGPz_kp1 = [x_k; x_kp1; nab_zero];                %inital vector for DPGPDA

DGEz_kp1 = [x_zero; -nab_zero; nab_zero];          %inital vector for DGT

DATCGTz_kp1 = [x_zero; -nab_zero; nab_zero];       %inital vector for DATCGT

DGTz_kp1 = [x_zero; -DGTeta1*nab_zero; nab_zero];   %inital vector for DGT

DooGTz_kp1 = [x_zero; -nab_zero; nab_zero];       %inital vector for DOO-GT

DGTgenz_kp1 = [x_zero; -(gamma+delta)*nab_zero; nab_zero];      %inital vector for oGGT

for ts=1:max_t
    ts
    %DEXTRA
    M_1=zeros(nodes,nodes);
    M_2=eye(nodes);
    M_3=-0.5*(eye(nodes)+W);
    M_4= 1*(eye(nodes)+W);
    M_5= DExeta1*eye(nodes);
    F=[-1*ones(1,nodes) 1*ones(1,nodes) DExeta1*ones(1,nodes)];
    G=zeros(1,nodes);
    %%% updating algorithm 
    DExz_k = DExz_kp1;
    [DExz_kp1,check] = GGT(DExz_k,M_1,M_2,M_3,M_4,M_5,F,G,nodes,2*trgt,y_data(:,ts+1),C_data);
    DEx_x_iter(ts,:) = DExz_kp1(1:6*nodes); %all node all signals
    


    %DNIDS
    M_1=zeros(nodes,nodes);
    M_2=eye(nodes);
    M_3=-0.5*(eye(nodes)+W);
    M_4= (eye(nodes)+W);
    M_5= DNseta1*0.5*(eye(nodes)+W);
    F=[-1*ones(1,nodes) 1*ones(1,nodes) DNseta1*ones(1,nodes)];
    G=zeros(1,nodes);
    %%% updating algorithms 
    DNsz_k = DNsz_kp1;
    [DNsz_kp1,check] = GGT(DNsz_k,M_1,M_2,M_3,M_4,M_5,F,G,nodes,2*trgt,y_data(:,ts+1),C_data);
    DNs_x_iter(ts,:) = DNsz_kp1(1:6*nodes);  %all node all signals



    %DPGPDA
    M_1=zeros(nodes,nodes);
    M_2=eye(nodes);
    M_3=-0.5*(eye(nodes)+Z_DPGP);
    M_4= eye(nodes)+Z_DPGP;
    M_5= 0.5*(DPGPc1)*inv(diag(d));
    F=[-1*ones(1,nodes) 1*ones(1,nodes) 0.5*(DPGPc1)*ones(1,nodes)];
    G=zeros(1,nodes);
    %%% updating algorithms 
    DPGPz_k = DPGPz_kp1;
    [DPGPz_kp1,check] = GGT(DPGPz_k,M_1,M_2,M_3,M_4,M_5,F,G,nodes,2*trgt,y_data(:,ts+1),C_data);
    DPGPDA_x_iter(ts,:) = DPGPz_kp1(1:6*nodes); %all node all signals


    
   % DGE
    M_1=W;
    M_2= DGEeta1*eye(nodes);
    M_3=zeros(nodes);
    M_4= W;
    M_5= eye(nodes);
    F=[zeros(1,nodes) 1*ones(1,nodes) ones(1,nodes)];
    G=zeros(1,nodes);
    %%% updating algorithms
    DGEz_k = DGEz_kp1;
    [DGEz_kp1,check] = GGT(DGEz_k,M_1,M_2,M_3,M_4,M_5,F,G,nodes,2*trgt,y_data(:,ts+1),C_data);
    DGE_x_iter(ts,:) = DGEz_kp1(1:6*nodes); %all node all signals


    %DATCGT
    M_1=W;
    M_2=DATCGTeta1 * W;
    M_3=zeros(nodes);
    M_4= W;
    M_5= W;
    F=[zeros(1,nodes) 1*ones(1,nodes) ones(1,nodes)];
    G=zeros(1,nodes);
    %%% updating algorithms
    DATCGTz_k = DATCGTz_kp1;
    [DATCGTz_kp1,check] = GGT(DATCGTz_k,M_1,M_2,M_3,M_4,M_5,F,G,nodes,2*trgt,y_data(:,ts+1),C_data);
    DATCGT_x_iter(ts,:) = DATCGTz_kp1(1:6*nodes); %all node all signals


   %DGT
    M_1=W;
    M_2=W;
    M_3=zeros(nodes);
    M_4= W;
    M_5= DGTeta1*eye(nodes);
    F=[zeros(1,nodes) 1*ones(1,nodes) DGTeta1*ones(1,nodes)];
    G=zeros(1,nodes);
    %%% updating algorithms
    DGTz_k = DGTz_kp1;
    [DGTz_kp1,check] = GGT(DGTz_k,M_1,M_2,M_3,M_4,M_5,F,G,nodes,2*trgt,y_data(:,ts+1),C_data);
    DGT_x_iter(ts,:) = DGTz_kp1(1:6*nodes); %all node all signals

   % DooGT
    M_1=W;
    M_2=DooGTeta1*W;
    M_3=zeros(nodes);
    M_4= W;
    M_5= eye(nodes);
    F=[zeros(1,nodes) 1*ones(1,nodes) 1*ones(1,nodes)];
    G=zeros(1,nodes);
    %%% updating algorithms
    DooGTz_k = DooGTz_kp1;
    [DooGTz_kp1,check] = GGT(DooGTz_k,M_1,M_2,M_3,M_4,M_5,F,G,nodes,2*trgt,y_data(:,ts+1),C_data);
    DooGT_x_iter(ts,:) = DooGTz_kp1(1:6*nodes); %all node all signals
    
    % opt-GGT
    M_1=W;
    M_2= alpha*eye(nodes)+ beta*W;
    M_3= zeros(nodes);
    M_4= W;
    M_5= gamma*eye(nodes) + delta*W;
    F=[zeros(1,nodes) 1*ones(1,nodes) ones(1,nodes)*M_5];
    G=zeros(1,nodes);
    %%% updating algorithms
    DGTgenz_k = DGTgenz_kp1;
    [DGTgenz_kp1,check] = GGT(DGTgenz_k,M_1,M_2,M_3,M_4,M_5,F,G,nodes,2*trgt,y_data(:,ts+1),C_data);
    DGTgen_x_iter(ts,:) = DGTgenz_kp1(1:6*nodes);  %all node all signals
end

%% finding Regret

%DEXTRA
fn_diff_DEx=zeros(max_t,1);
for j=1:nodes
    for k=1:max_t
        fn_diff_DEx(k) = fn_diff_DEx(k) + 0.5*norm(C_data*DEx_x_iter(k,(j-1)*2*trgt+1:j*2*trgt)'-y_data(:,k+1));
    end
end
fn_diff_DEx = (1/nodes)*fn_diff_DEx;

for k=1:max_t
    f_star(k) = 0.5*norm(C_data*x_data(:,k+1)-y_data(:,k+1));
end

Reg_K_DEx(1)=fn_diff_DEx(1);
for i=2:max_t
    Reg_K_DEx(i)=fn_diff_DEx(i)+Reg_K_DEx(i-1);
end



%DNIDS
fn_diff_DNs=zeros(max_t,1);
for j=1:nodes
    for k=1:max_t
        fn_diff_DNs(k) = fn_diff_DNs(k) + 0.5*norm(C_data*DNs_x_iter(k,(j-1)*2*trgt+1:j*2*trgt)'-y_data(:,k+1));
    end
end
fn_diff_DNs = (1/nodes)*fn_diff_DNs;

Reg_K_DNs(1)=fn_diff_DNs(1);
for i=2:max_t
    Reg_K_DNs(i)=fn_diff_DNs(i)+Reg_K_DNs(i-1);
end


%DPGPDA
fn_diff_DPGPDA=zeros(max_t,1);
for j=1:nodes
    for k=1:max_t
        fn_diff_DPGPDA(k) = fn_diff_DPGPDA(k) + 0.5*norm(C_data*DPGPDA_x_iter(k,(j-1)*2*trgt+1:j*2*trgt)'-y_data(:,k+1));
    end
end
fn_diff_DPGPDA = (1/nodes)*fn_diff_DPGPDA;

Reg_K_DPGPDA(1)=fn_diff_DPGPDA(1);
for i=2:max_t
    Reg_K_DPGPDA(i)=fn_diff_DPGPDA(i)+Reg_K_DPGPDA(i-1);
end

%DGE
fn_diff_DGE=zeros(max_t,1);
for j=1:nodes
    for k=1:max_t
        fn_diff_DGE(k) = fn_diff_DGE(k) + 0.5*norm(C_data*DGE_x_iter(k,(j-1)*2*trgt+1:j*2*trgt)'-y_data(:,k+1));
    end
end
fn_diff_DGE = (1/nodes)*fn_diff_DGE;

Reg_K_DGE(1)=fn_diff_DGE(1);
for i=2:max_t
    Reg_K_DGE(i)=fn_diff_DGE(i)+Reg_K_DGE(i-1);
end

%DATCGT
fn_diff_DATCGT=zeros(max_t,1);
for j=1:nodes
    for k=1:max_t
        fn_diff_DATCGT(k) = fn_diff_DATCGT(k) + 0.5*norm(C_data*DATCGT_x_iter(k,(j-1)*2*trgt+1:j*2*trgt)'-y_data(:,k+1));
    end
end
fn_diff_DATCGT = (1/nodes)*fn_diff_DATCGT;

Reg_K_DATCGT(1)=fn_diff_DATCGT(1);
for i=2:max_t
    Reg_K_DATCGT(i)=fn_diff_DATCGT(i)+Reg_K_DATCGT(i-1);
end

%DGT
fn_diff_DGT=zeros(max_t,1);
for j=1:nodes
    for k=1:max_t
        fn_diff_DGT(k) = fn_diff_DGT(k) + 0.5*norm(C_data*DGT_x_iter(k,(j-1)*2*trgt+1:j*2*trgt)'-y_data(:,k+1));
    end
end
fn_diff_DGT = (1/nodes)*fn_diff_DGT;

Reg_K_DGT(1)=fn_diff_DGT(1);
for i=2:max_t
    Reg_K_DGT(i)=fn_diff_DGT(i)+Reg_K_DGT(i-1);
end

%DooGT
fn_diff_DooGT=zeros(max_t,1);
for j=1:nodes
    for k=1:max_t
        fn_diff_DooGT(k) = fn_diff_DooGT(k) + 0.5*norm(C_data*DooGT_x_iter(k,(j-1)*2*trgt+1:j*2*trgt)'-y_data(:,k+1));
    end
end
fn_diff_DooGT = (1/nodes)*fn_diff_DooGT;

Reg_K_DooGT(1)=fn_diff_DooGT(1);
for i=2:max_t
    Reg_K_DooGT(i)=fn_diff_DooGT(i)+Reg_K_DooGT(i-1);
end

%oGGT
fn_diff_DGTgen=zeros(max_t,1);
for j=1:nodes
    for k=1:max_t
        fn_diff_DGTgen(k) = fn_diff_DGTgen(k) + 0.5*norm(C_data*DGTgen_x_iter(k,(j-1)*2*trgt+1:j*2*trgt)'-y_data(:,k+1));
    end
end
fn_diff_DGTgen = (1/nodes)*fn_diff_DGTgen;

Reg_K_DGTgen(1)=fn_diff_DGTgen(1);
for i=2:max_t
    Reg_K_DGTgen(i)=fn_diff_DGTgen(i)+Reg_K_DGTgen(i-1);
end


%% Plotting Results 
% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#D95319';
color1 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

figure(2)
plot(ratio,DExrho,"Color",'r','linewidth',2);
xlabel('stepsize(\eta)*L','FontSize',12,'FontWeight','bold')
ylabel('\rho','FontSize',12,'FontWeight','bold')
hold on
plot(ratio,DNsrho,"Color",'g','linewidth',2);
plot(ratio,DPGPrho,"Color",'c','linewidth',2);
plot(ratio,DGErho,"Color",'m','linewidth',2);
plot(ratio,DATrho,"Color",'y','linewidth',2);
plot(ratio,DGTrho,"Color",'b','linewidth',2);
plot(ratio,DooGTrho,"Color",color1,'linewidth',2);
plot(ratio,oGGTrho,'--',"Color",'k','linewidth',2);
legend('DEXTRA','DNIDS','DPGPDA','DGE','DATCGT','DGT','DOOGT','oGGT','FontSize',10,'FontWeight','bold');
grid on




% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#D95319';
color1 = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

figure(3)
semilogy(1:max_t, Reg_K_DEx,'linewidth',2,"color",'r')
grid on
xlabel('K (iterations)','FontSize',14,'FontWeight','bold')
ylabel('Reg^K','FontSize',14,'FontWeight','bold')
hold on 
semilogy(1:max_t, Reg_K_DNs,'linewidth',2,"color",'g')
semilogy(1:max_t, Reg_K_DPGPDA,'linewidth',2,"color",'c')
semilogy(1:max_t, Reg_K_DGE,'linewidth',2,"color",'m')
semilogy(1:max_t, Reg_K_DATCGT,'linewidth',2,"color",'y')
semilogy(1:max_t, Reg_K_DGT,'linewidth',2,"color",'b')
semilogy(1:max_t, Reg_K_DooGT,'linewidth',2,"color",color1)
semilogy(1:max_t, Reg_K_DGTgen,'linewidth',2,"color",'k')
legend('DEXTRA','DNIDS','DPGPDA','DGE','DATCGT','DGT','DOOGT','oGGT','FontSize',10,'FontWeight','bold');

