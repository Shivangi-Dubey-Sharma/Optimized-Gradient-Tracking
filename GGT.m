function [z_kp1,check] = GGT(z_k,M_1,M_2,M_3,M_4,M_5,F,G,n,d,y_data,C_data)
%This function gives the next iterate value of GGT algorithm for considered
%problem
%   Detailed explanation goes here
% inputs
% z_k = [x_km1, s_k, nab_km1]                   vector  (3nd X 1)
% M_1                                           matrix  (n X n)
% M_2                                           matrix  (n X n)
% M_3                                           matrix  (n X n)
% M_4                                           matrix  (n X n)
% M_5                                           matrix  (n X n)
% F                                             vector  (1 X 3n)
% n = number of nodes in the network
% d = dimension of individual node vector x_k   scalar  
% y_data = observation
% C_data = measurement matrix, C^TC>0           matrix  (d X d)
%outputs
% z_kp1 = next iterate vector (state vector)    vector  (3nd X 1)
% check = just to verify the iterate vector(22(c))     vector  zeros(d X 1)

%% defining A B C and D
A=[M_1          M_2      zeros(n,n);
   M_3          M_4         M_5;
   zeros(n,n) zeros(n,n) zeros(n,n) ];

B=[zeros(n,n);
    -M_5;
    eye(n,n)];

C=[M_1 M_2 zeros(n,n)];

D=zeros(n,n);


y_k = kron(C,eye(d))*z_k;
nab = [];
for i=1:n
    i;
    nab = [nab; C_data((i-1)*6+1:i*6,:)'*(C_data((i-1)*6+1:i*6,:)*y_k((i-1)*d+1:i*d,1) - y_data((i-1)*6+1:i*6,1))];
end
u_k = nab;

z_kp1 = kron(A,eye(d))*z_k + kron(B,eye(d))*u_k; 
check = kron(F,eye(d))*z_k ; 
end