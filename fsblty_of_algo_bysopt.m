function [outres] = fsblty_of_algo_bysopt(n,d,Mfnl,M_1,M_2,M_3,M_4,M_5,F,G,ro_vec)
% fsblty_of_algo function provides the matrix P and rho for which the
% following SDP is feasible 
% [A'PA - rho^2 P    A'PB            [C  D      [C  D
%  B'PA              B'PB ] + lambda  0  I]' M   0  I] <= 0
% where M is 
% M = [diag({-2 L_i*(L_i*I_d - M_i)^(-1)*M_i}_{i=1:n})          diag({*(L_i*I_d - M_i)^(-1)*(L_i*I_d + M_i)}_{i=1:n})
%       diag({*(L_i*I_d - M_i)^(-1)*(L_i*I_d + M_i)}_{i=1:n})         -2diag({*(L_i*I_d - M_i)^(-1)}_{i=1:n})           ] \in R^(2nd X 2nd)
%   Detailed explanation goes here

%% inputs
% n : Number of nodes in the network
% d : Dimention of optimization variable
% Mfnl : Matrix M of SDP in paper ( equation 25) 
% M_1 to M_5 are matrices defining special case of GGT \in R^(n X n)
% F matrices defining special case of GGT \in R^(1 X 3n)
% G matrices defining special case of GGT \in R^(1 X n)

%% Outputs
%Pf : P matrix of SDP \in R^(3nd X 3nd) for minimum rho of given algorithm
%rho : minimum \rho of the considered special case of GGT
%lambda1 : lambda of SDP for above \rho and considered algorithm
%Cvx_s : Status of CVX for all \rho evaluated in the bisection
%t : Condition number of Matrix P


%%
M=Mfnl;

%% defining A B C and D of SDP
A=[M_1          M_2      zeros(n,n);
   M_3          M_4         M_5;
   zeros(n,n) zeros(n,n) zeros(n,n) ];
Af = kron(A,eye(d));


B=[zeros(n,n);
    -M_5;
    eye(n,n)];
Bf = kron(B,eye(d));


C=[M_1 M_2 zeros(n,n)];
Cf = kron(C,eye(d));

D=zeros(n,n);
Df = kron(D,eye(d));

R=null([F G]);
Rf = null(kron([F G],eye(d)));


Cvx_s=strings;

%% feasible SDP parameter finding %%

%% bisection
% bisection is done to find minimum \rho 
% rol : lower limit of the range of rho
% rou : upper limit of the range of rho

rol = ro_vec(1); 
rou = ro_vec(2);
num_run = 0;
checki = 0;     %to keep running the bisection algorithm even though rou-rol is small,
                % until we get a feasible point
ro = rol+0.5*(rou-rol);  % starting point
solvi = 0;      % solvi = 1 represents that we get a feasible point

%muu = 5*10^-5;   % As CVX does not take Positive Definite Matrix Condition we 
                % have made the condition of P to be PSD and added muu*I with P, 
                % wherever P is defined in SDP 

%muu = 5*10^-6;

while (rou-rol>10^-7 || checki==0)
 
    num_run=num_run+1;  %just to see the current position of simulation
    ro ;                %to see for which ro value the code is running     


% updated form of SDP to minimize the condition number of P 
     cvx_begin sdp quiet
        variable P(3*n*d,3*n*d) symmetric 
        variables lmbdatelda t
        minimize( t );
        subject to
        P <= t*eye(3*n*d,3*n*d);
        P >= eye(3*n*d,3*n*d);
        X = Rf'*([Af'*P*Af-ro^2*P Af'*P*Bf; Bf'*P*Af Bf'*P*Bf]+ lmbdatelda*[Cf Df; zeros(n*d,3*n*d) eye(n*d)]'*M*[Cf Df; zeros(n*d,3*n*d) eye(n*d,n*d)])*Rf;
        X+X' <= 0
        lmbdatelda>=0
    cvx_end


    %for strongly convex function kronecker product with d can be factored
    %out so use following code with simplified Mfnl ( R^(2n X 2n)) as input
%     cvx_begin sdp quiet
%         variable P(3*n,3*n) symmetric 
%         variables lmbdatelda t
%         minimize( t );
%         subject to
%         P <= t*eye(3*n,3*n);
%         P >= eye(3*n,3*n);
%         X = R'*([A'*P*A-ro^2*P A'*P*B; B'*P*A B'*P*B]+ lmbdatelda*[C D; zeros(n,3*n) eye(n)]'*M*[C D; zeros(n,3*n) eye(n,n)])*R;
%         X+X' <= 0
%         lmbdatelda>=0
%     cvx_end


    s_nan = sum(sum(isnan(P)))>0; % 1 if matrix P contains NAN 
    s_inf = sum(sum(isinf(P)))>0; % 1 if matrix P contains Inf

    if ((strcmp(cvx_status,'Solved')) && (lmbdatelda~=NaN) && (~s_nan) && (~s_inf))  %verifying the condition that we get a feasible point
       rou = ro;            % solved so search in the lower half
       checki = 1;          %to stop if its the last round (as rou-rol is very small)
       solvi = 1            %to check if its getting solved or not
       rho = ro;            %updating Final rho
       lmbda1 = lmbdatelda;   %updating Final lambda
       Cvx_s = [Cvx_s,cvx_status];         % Storing CVX staus 
       Pf = P;
    else 
        rol = ro;      % not solved so search in the upper half
        checki = 0;    %to move on if its still infeasible even if its a last round
    end
    ro = rol + 0.5*(rou-rol);  % updating the mid point

    if num_run>500           % to stop it from going in an infinite loop
        if solvi==0         % if problem is still infeasible
            rho = 1.1;
            lmbda1= lmbdatelda; %value of last round (we can write any random value)
            Cvx_s = [Cvx_s,cvx_status];
            Pf = P;
        end
        break
    end
end 
outres.factrfnl = ((rho.*rho)./((1-rho).^2)).*t;
outres.Pf = Pf;
outres.rho = rho;
outres.lmbda1 = lmbda1; 
outres.status = Cvx_s;
outres.t = t;
    
 end

