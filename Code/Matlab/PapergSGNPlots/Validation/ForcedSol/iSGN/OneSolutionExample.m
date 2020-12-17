% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolSolitonLoop/06/";

wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/ForcedSolutions/iSGN/05/"

param = fileread(strcat(wdir ,'Param.dat' ));
dxstr = extractBetween(param,"(x_len -1)  :"," n_GhstCells");
dx = str2double(dxstr{1,1});

scenbeg = -100;
scenend = 100;

Region1LB = -25;
Region1UB = 38;
Region2UB = 62;
Region3UB = 75;


Region1LB0 = -25;
Region1UB0 = -15;
Region2UB0 = 15;
Region3UB0 = 75;

denselinesep = 1;
densedotsep = 10;
lessdenselinesep = 10;
lessdensedotsep = 80;

Region1LBi = round(((Region1LB - scenbeg)/dx));
Region1UBi = round(((Region1UB - scenbeg)/dx));
Region2UBi  = round(((Region2UB - scenbeg)/dx));
Region3UBi = round(((Region3UB  - scenbeg)/dx));


Region1LB0i = round(((Region1LB0 - scenbeg)/dx));
Region1UB0i = round(((Region1UB0 - scenbeg)/dx));
Region2UB0i  = round(((Region2UB0 - scenbeg)/dx));
Region3UB0i = round(((Region3UB0  - scenbeg)/dx));



End = importdata(strcat(wdir,'End.dat' ));
Init = importdata(strcat(wdir,'Init.dat' ));

t0 = Init(1,1);
x0 = [Init(Region1LB0i:lessdenselinesep:Region1UB0i,2); ...
    Init(Region1UB0i:denselinesep:Region2UB0i,2); ...
    Init(Region2UB0i:lessdenselinesep:Region3UB0i,2)];
h0 = [Init(Region1LB0i:lessdenselinesep:Region1UB0i,3); ...
    Init(Region1UB0i:denselinesep:Region2UB0i,3); ...
    Init(Region2UB0i:lessdenselinesep:Region3UB0i,3)];
G0 = [Init(Region1LB0i:lessdenselinesep:Region1UB0i,4); ...
    Init(Region1UB0i:denselinesep:Region2UB0i,4); ...
    Init(Region2UB0i:lessdenselinesep:Region3UB0i,4)];
u0 = [Init(Region1LB0i:lessdenselinesep:Region1UB0i,5); ...
    Init(Region1UB0i:denselinesep:Region2UB0i,5); ...
    Init(Region2UB0i:lessdenselinesep:Region3UB0i,5)];


t = End(1,1);
x = [End(Region1LBi:lessdensedotsep:Region1UBi,2); ...
    End(Region1UBi:densedotsep:Region2UBi,2); ...
    End(Region2UBi:lessdensedotsep:Region3UBi,2)];
h = [End(Region1LBi:lessdensedotsep:Region1UBi,3); ...
    End(Region1UBi:densedotsep:Region2UBi,3); ...
    End(Region2UBi:lessdensedotsep:Region3UBi,3)];
G = [End(Region1LBi:lessdensedotsep:Region1UBi,5); ...
    End(Region1UBi:densedotsep:Region2UBi,5); ...
    End(Region2UBi:lessdensedotsep:Region3UBi,5)];
u = [End(Region1LBi:lessdensedotsep:Region1UBi,7); ...
    End(Region1UBi:densedotsep:Region2UBi,7); ...
    End(Region2UBi:lessdensedotsep:Region3UBi,7)];


tA = End(1,1);
xA = [End(Region1LBi:lessdenselinesep:Region1UBi,2); ...
    End(Region1UBi:denselinesep:Region2UBi,2); ...
    End(Region2UBi:lessdenselinesep:Region3UBi,2)];
hA = [End(Region1LBi:lessdenselinesep:Region1UBi,4); ...
    End(Region1UBi:denselinesep:Region2UBi,4); ...
    End(Region2UBi:lessdenselinesep:Region3UBi,4)];
GA = [End(Region1LBi:lessdenselinesep:Region1UBi,6); ...
    End(Region1UBi:denselinesep:Region2UBi,6); ...
    End(Region2UBi:lessdenselinesep:Region3UBi,6)];
uA = [End(Region1LBi:lessdenselinesep:Region1UBi,8); ...
    End(Region1UBi:denselinesep:Region2UBi,8); ...
    End(Region2UBi:lessdenselinesep:Region3UBi,8)];



figure;
plot(x,h,'.r');
hold on;
plot(xA,hA,'-b');
plot(x0,h0,'-k');
hold off
xlabel('x (m)');
ylabel('h (m)');
axis([-25 75 0.9,1.6]);
xticks([-25,0,25,50,75]);
yticks([0.9,1,1.1,1.2,1.3,1.4,1.5,1.6]);
legend('hide')
matlab2tikz('h.tex');
close all;


figure;
plot(x,G,'.r');
hold on;
plot(xA,GA,'-b');
plot(x0,G0,'-k');
hold off
xlabel('x (m)');
ylabel('G (m^2/s)');
axis([-25 75 -0.1 0.5]);
xticks([-25,0,25,50,75]);
yticks([-0.1,0,0.1,0.2,0.3,0.4,0.5]);
legend('hide')
matlab2tikz('G.tex');
close all;

figure;
plot(x,u,'.r');
hold on;
plot(xA,uA,'-b');
plot(x0,u0,'-k');
hold off
xlabel('x (m)');
% ylabel('u (m/s)');
axis([-25 75 -0.1 0.4]);
xticks([-25,0,25,50,75]);
yticks([-0.1,0,0.1,0.2,0.3,0.4]);
legend('hide');
matlab2tikz('u.tex');


