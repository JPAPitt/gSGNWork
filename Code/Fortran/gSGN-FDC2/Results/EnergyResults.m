% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data

 wdir = "./FD2C/Validation/AnalyticSolutions/SerreSoliton/"



OutEps = importdata(strcat(wdir, 'Energy.dat'));
dxseps = OutEps(:,1);
Ch = OutEps(:,2);
Cuh = OutEps(:,3);
CE = OutEps(:,4);

figure;
loglog(dxseps,Ch,'s b',dxseps,Cuh,'^ k',dxseps,CE,'x g',dxseps,dxseps.^2,'--k')

grid off
legend('hide')
xlabel('\Delta x')
ylabel('L_2')




