% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/AnalyticSolutions/DBSWWE/"

OutEps = importdata(strcat(wdir, 'Energy.dat'));
dxs = OutEps(:,1);
Ehs = OutEps(:,2);
EGs = OutEps(:,3);
Euhs = OutEps(:,4);
EHs = OutEps(:,5);

figure;
loglog(dxs,Ehs,'s b',dxs,EGs,'o r',dxs,Euhs,'^ k',dxs,EHs,'x g')
grid off
legend('h','G', 'uh', 'H','Location','northwest')
xlabel('\Delta x')
ylabel('C_1')
axis([10^-4 10 10^-16 1]);
xticks([10^-4,10^-3,10^-2,10^-1,10^0,10]);
yticks([10^-16,10^-12,10^-8,10^-4,1]);
% matlab2tikz('EnergyResults.tex');




