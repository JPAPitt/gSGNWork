% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/AnalyticSolutions/DBSWWETheta1/05/"


param = fileread(strcat(wdir ,'Params.dat' ));
dxstr = extractBetween(param,"dx","tstart");
dx = str2double(dxstr{1,1});

scenbeg = -250;
scenend = 250;


Region1LB = -175;
Region1UB = -155;
Region2UB = -80;
Region3UB = 147;
Region4UB = 151;
Region5UB = 175;

denselinesep = 1;
densedotsep = 1;
middenselinesep = 1;
middensedotsep = 1;
lessdenselinesep = 1;
lessdensedotsep = 1;



Region1LBi = round(((Region1LB - scenbeg)/dx));
Region1UBi = round(((Region1UB - scenbeg)/dx));
Region2UBi  = round(((Region2UB - scenbeg)/dx));
Region3UBi = round(((Region3UB  - scenbeg)/dx));
Region4UBi= round(((Region4UB  - scenbeg)/dx));
Region5UBi= round(((Region5UB  - scenbeg)/dx));

EndA = importdata(strcat(wdir,'EndAna.dat' ));
End = importdata(strcat(wdir,'End.dat' ));
Init = importdata(strcat(wdir,'Init.dat' ));


t0 = Init(1,1);
x0 = [Init(Region1LBi:lessdenselinesep:Region1UBi,2); ...
    Init(Region1UBi:middenselinesep:Region2UBi,2); ...
    Init(Region2UBi:lessdenselinesep:Region3UBi,2); ...
    Init(Region3UBi:denselinesep:Region4UBi,2); ...
    Init(Region4UBi:lessdenselinesep:Region5UBi,2)];
h0 = [Init(Region1LBi:lessdenselinesep:Region1UBi,3); ...
    Init(Region1UBi:middenselinesep:Region2UBi,3); ...
    Init(Region2UBi:lessdenselinesep:Region3UBi,3); ...
    Init(Region3UBi:denselinesep:Region4UBi,3); ...
    Init(Region4UBi:lessdenselinesep:Region5UBi,3)];
G0 = [Init(Region1LBi:lessdenselinesep:Region1UBi,4); ...
    Init(Region1UBi:middenselinesep:Region2UBi,4); ...
    Init(Region2UBi:lessdenselinesep:Region3UBi,4); ...
    Init(Region3UBi:denselinesep:Region4UBi,4); ...
    Init(Region4UBi:lessdenselinesep:Region5UBi,4)];
u0 = [Init(Region1LBi:lessdenselinesep:Region1UBi,5); ...
    Init(Region1UBi:middenselinesep:Region2UBi,5); ...
    Init(Region2UBi:lessdenselinesep:Region3UBi,5); ...
    Init(Region3UBi:denselinesep:Region4UBi,5); ...
    Init(Region4UBi:lessdenselinesep:Region5UBi,5)];


t = End(1,1);
x = [End(Region1LBi:lessdensedotsep:Region1UBi,2); ...
    End(Region1UBi:middensedotsep:Region2UBi,2); ...
    End(Region2UBi:lessdensedotsep:Region3UBi,2); ...
    End(Region3UBi:densedotsep:Region4UBi,2); ...
    End(Region4UBi:lessdensedotsep:Region5UBi,2)];
h = [End(Region1LBi:lessdensedotsep:Region1UBi,3); ...
    End(Region1UBi:middensedotsep:Region2UBi,3); ...
    End(Region2UBi:lessdensedotsep:Region3UBi,3); ...
    End(Region3UBi:densedotsep:Region4UBi,3); ...
    End(Region4UBi:lessdensedotsep:Region5UBi,3)];
G = [End(Region1LBi:lessdensedotsep:Region1UBi,4); ...
    End(Region1UBi:middensedotsep:Region2UBi,4); ...
    End(Region2UBi:lessdensedotsep:Region3UBi,4); ...
    End(Region3UBi:densedotsep:Region4UBi,4); ...
    End(Region4UBi:lessdensedotsep:Region5UBi,4)];
u = [End(Region1LBi:lessdensedotsep:Region1UBi,5); ...
    End(Region1UBi:middensedotsep:Region2UBi,5); ...
    End(Region2UBi:lessdensedotsep:Region3UBi,5); ...
    End(Region3UBi:densedotsep:Region4UBi,5); ...
    End(Region4UBi:lessdensedotsep:Region5UBi,5)];


tA = EndA(1,1);
xA = [EndA(Region1LBi:lessdenselinesep:Region1UBi,2); ...
    EndA(Region1UBi:middenselinesep:Region2UBi,2); ...
    EndA(Region2UBi:lessdenselinesep:Region3UBi,2); ...
    EndA(Region3UBi:denselinesep:Region4UBi,2); ...
    EndA(Region4UBi:lessdenselinesep:Region5UBi,2)];
hA = [EndA(Region1LBi:lessdenselinesep:Region1UBi,3); ...
    EndA(Region1UBi:middenselinesep:Region2UBi,3); ...
    EndA(Region2UBi:lessdenselinesep:Region3UBi,3); ...
    EndA(Region3UBi:denselinesep:Region4UBi,3); ...
    EndA(Region4UBi:lessdenselinesep:Region5UBi,3)];
GA = [EndA(Region1LBi:lessdenselinesep:Region1UBi,4); ...
    EndA(Region1UBi:middenselinesep:Region2UBi,4); ...
    EndA(Region2UBi:lessdenselinesep:Region3UBi,4); ...
    EndA(Region3UBi:denselinesep:Region4UBi,4); ...
    EndA(Region4UBi:lessdenselinesep:Region5UBi,4)];
uA = [EndA(Region1LBi:lessdenselinesep:Region1UBi,5); ...
    EndA(Region1UBi:middenselinesep:Region2UBi,5); ...
    EndA(Region2UBi:lessdenselinesep:Region3UBi,5); ...
    EndA(Region3UBi:denselinesep:Region4UBi,5); ...
    EndA(Region4UBi:lessdenselinesep:Region5UBi,5)];



figure;
plot(x,h,'.r');
hold on;
plot(xA,hA,'-b');
plot(x0,h0,'-k');
hold off
xlabel('x (m)');
ylabel('h (m)');
axis([-175 175  0.9 2.1]);
xticks([-175,-87.5,0,87.5,175]);
yticks([0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1]);
legend('hide')
% matlab2tikz('hEx.tex');
% close all;


% figure;
% plot(x,G,'.r');
% hold on;
% plot(xA,GA,'-b');
% plot(x0,G0,'-k');
% hold off
% xlabel('x (m)');
% ylabel('G (m^2/s)');
% axis([-175 175  0 2]);
% xticks([-175,-87.5,0,87.5,175]);
% yticks([0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2]);
% legend('hide')
% matlab2tikz('GEx.tex');
% close all;
% 
% 
% figure;
% plot(x,u,'.r');
% hold on;
% plot(xA,uA,'-b');
% plot(x0,u0,'-k');
% hold off
% xlabel('x (m)');
% ylabel('u (m/s)');
% axis([-175 175  0 1.4]);
% xticks([-175,-87.5,0,87.5,175]);
% yticks([0,0.2,0.4,0.6,0.8,1,1.2,1.4]);
% legend('hide');
% matlab2tikz('uEx.tex');
% close all;

