function [W,G,d] = MH_algo_comm_mat(p,n)
%DISTRIBUTED METROPOLIS–HASTINGS ALGORITHM
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9239886&tag=1
% input agr : p is the probability of having the connection between nodes, 
% small p implies sparse network 
% n is the number of nodes in the network.
% output arg : W is the doubly stochastic symmetric matrix
% output arg : G is a symmetric matrix representing Graph (1 on connected node, 0 on diagonals)
% output arg : d is the degree of each node (column vector)
for i=1:n
    for j=1:i
        if i==j
            G(i,j)=0;
        else
            G(i,j)=binornd(1,p);
            G(j,i) = G(i,j);
        end
    end
end

for i=1:n
    d(i)=sum(G(i,:)); %degree of each node i
end

for i=1:n
    N(i).node =  find(G(i,:)); %%%% structure n stores the name of nebouring node in field node
end

W=zeros(n,n);
[roG1,colG1] = find(G);
for i=1:length(roG1)
    W(roG1(i),colG1(i)) = 1/max(d(roG1(i)),d(colG1(i)));
end
for i=1:n
    W(i,i) = 1-sum(W(i,:));
end

end