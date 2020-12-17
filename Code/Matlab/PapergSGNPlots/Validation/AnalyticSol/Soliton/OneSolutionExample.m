% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolSolitonLoop/06/";
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/AnalyticSolutions/SerreSoliton/06/"


param = fileread(strcat(wdir ,'Params.dat' ));
dxstr = extractBetween(param,"dx","tstart");
dx = str2double(dxstr{1,1});

scenbeg = -200;
scenend = 200;


Region1LB = -25;
Region1UB = 118;
Region2UB = 128;
Region3UB = 150;


Region1LB0 = -25;
Region1UB0 = -10;
Region2UB0 = 10;
Region3UB0 = 150;

denselinesep = 1;
densedotsep = 6;
lessdenselinesep = 10;
lessdensedotsep = 60;

Region1LBi = round(((Region1LB - scenbeg)/dx));
Region1UBi = round(((Region1UB - scenbeg)/dx));
Region2UBi  = round(((Region2UB - scenbeg)/dx));
Region3UBi = round(((Region3UB  - scenbeg)/dx));


Region1LB0i = round(((Region1LB0 - scenbeg)/dx));
Region1UB0i = round(((Region1UB0 - scenbeg)/dx));
Region2UB0i  = round(((Region2UB0 - scenbeg)/dx));
Region3UB0i = round(((Region3UB0  - scenbeg)/dx));



EndA = importdata(strcat(wdir,'EndAna.dat' ));
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
G = [End(Region1LBi:lessdensedotsep:Region1UBi,4); ...
    End(Region1UBi:densedotsep:Region2UBi,4); ...
    End(Region2UBi:lessdensedotsep:Region3UBi,4)];
u = [End(Region1LBi:lessdensedotsep:Region1UBi,5); ...
    End(Region1UBi:densedotsep:Region2UBi,5); ...
    End(Region2UBi:lessdensedotsep:Region3UBi,5)];


tA = EndA(1,1);
xA = [EndA(Region1LBi:lessdenselinesep:Region1UBi,2); ...
    EndA(Region1UBi:denselinesep:Region2UBi,2); ...
    EndA(Region2UBi:lessdenselinesep:Region3UBi,2)];
hA = [EndA(Region1LBi:lessdenselinesep:Region1UBi,3); ...
    EndA(Region1UBi:denselinesep:Region2UBi,3); ...
    EndA(Region2UBi:lessdenselinesep:Region3UBi,3)];
GA = [EndA(Region1LBi:lessdenselinesep:Region1UBi,4); ...
    EndA(Region1UBi:denselinesep:Region2UBi,4); ...
    EndA(Region2UBi:lessdenselinesep:Region3UBi,4)];
uA = [EndA(Region1LBi:lessdenselinesep:Region1UBi,5); ...
    EndA(Region1UBi:denselinesep:Region2UBi,5); ...
    EndA(Region2UBi:lessdenselinesep:Region3UBi,5)];

figure;
plot(x,h,'.r');
hold on;
plot(xA,hA,'-b');
plot(x0,h0,'-k');
hold off
xlabel('x (m)');
ylabel('h (m)');
axis([-25 150 0.9,1.8]);
xticks([-25,0,25,50,75,100,125,150]);
yticks([0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8]);
legend('hide')
matlab2tikz('hEx.tex');
close all;


figure;
plot(x,G,'.r');
hold on;
plot(xA,GA,'-b');
plot(x0,G0,'-k');
hold off
xlabel('x (m)');
ylabel('G (m^2/s)');
axis([-25 150 -0.5 4]);
xticks([-25,0,25,50,75,100,125,150]);
yticks([-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4]);
legend('hide')
matlab2tikz('GEx.tex');
close all;

figure;
plot(x,u,'.r');
hold on;
plot(xA,uA,'-b');
plot(x0,u0,'-k');
hold off
xlabel('x (m)');
ylabel('u (m/s)');
axis([-25 150 -0.2 1.8]);
xticks([-25,0,25,50,75,100,125,150]);
yticks([-0.2,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8]);
legend('hide');
matlab2tikz('uEx.tex');
close all;

