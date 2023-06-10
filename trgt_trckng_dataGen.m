function [x_star,y_t] = trgt_trckng_dataGen(nodes,targets,amp_min,amp_max,periodicity,frequency,timestamp,C)
%UNTITLED3 Summary of this function goes here
%this function generates the data of target tracking example present in the following
%paper
%   Detailed explanation goes here
% seed = seed required for the genration of random variables    scalar
% nodes(n) = number of nodes in the network                     scalar
% targets(tr) = number of moving targets                        scalar
% amp_min = min value of the amplitude of targets               vector (1Xtr)
% amp_max = max value of amplitude                              vector (1Xtr)
% periodicity = periodicity of the signals to be tracked(sec)   vector (1Xtr)  
% frequency = sampling frequency  (Hz)                          vector (1Xtr)
% timestamp = which sample number                               scalar
% C= measurement matrix, C^TC>0                                 vector (2.trX2.tr)
%just keep in mind that p./f = n times ones(1 X tr) 



x_star=[];
for i=1:targets
    %A and phi are constant withrespect to time
    A =  amp_min(1,i) + (amp_max(1,i)-amp_min(1,i)).*rand(1,1);
    phi = pi.*rand(1,1);
    x_star = [x_star; A*sin(2*pi*(1/(frequency(1,i)* periodicity(1,i)))*timestamp + phi*ones(1,length(timestamp))); 
        2*pi/periodicity(1,i)*A*cos(2*pi*(1/(frequency(1,i)* periodicity(1,i)))*timestamp + phi*ones(1,length(timestamp)))];
end

y_t = C*x_star;

end