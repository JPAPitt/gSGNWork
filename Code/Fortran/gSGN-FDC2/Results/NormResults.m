% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data

 wdir = "./FD2C/Validation/AnalyticSolutions/SerreSoliton/"



OutEps = importdata(strcat(wdir, 'Norms.dat'));
dxseps = OutEps(:,1);
Normhseps = OutEps(:,2);
Normuseps = OutEps(:,3);

figure;
loglog(dxseps,Normhseps,'s b',dxseps,Normuseps,'^ k',dxseps,dxseps.^2,'--k')

grid off
legend('hide')
xlabel('\Delta x')
ylabel('L_2')




