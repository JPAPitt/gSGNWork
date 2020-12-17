% Process Fortran Outputs

clc;
clear all;
close all;

wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/Examples/DB/00/";

param = fileread(strcat(wdir ,'Params.dat' ));
dxstr = extractBetween(param,"dx","tstart");
dx = str2double(dxstr{1,1});

scenbeg = -250;
scenend = 250;


Region1LB00 = -200;
Region1UB00 = -100;
Region2UB00 = 150;
Region3UB00 = 200;


denselinesep = 5;
lessdenselinesep = 10;

Region1LBi00 = round(((Region1LB00 - scenbeg)/dx)) + 1;
Region1UBi00 = round(((Region1UB00 - scenbeg)/dx)) + 1;
Region2UBi00  = round(((Region2UB00 - scenbeg)/dx)) + 1;
Region3UBi00 = round(((Region3UB00  - scenbeg)/dx)) + 1;

End = importdata(strcat(wdir,'End.dat' ));


t = End(1,1);
x00 = [End(Region1LBi00:lessdenselinesep:Region1UBi00,2); ...
    End(Region1UBi00:denselinesep:Region2UBi00,2); ...
    End(Region2UBi00:lessdenselinesep:Region3UBi00,2)];
h00 = [End(Region1LBi00:lessdenselinesep:Region1UBi00,3); ...
    End(Region1UBi00:denselinesep:Region2UBi00,3); ...
    End(Region2UBi00:lessdenselinesep:Region3UBi00,3)];
G00 = [End(Region1LBi00:lessdenselinesep:Region1UBi00,4); ...
    End(Region1UBi00:denselinesep:Region2UBi00,4); ...
    End(Region2UBi00:lessdenselinesep:Region3UBi00,4)];
u00 = [End(Region1LBi00:lessdenselinesep:Region1UBi00,5); ...
    End(Region1UBi00:denselinesep:Region2UBi00,5); ...
    End(Region2UBi00:lessdenselinesep:Region3UBi00,5)];

wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/Examples/DB/01/";
End = importdata(strcat(wdir,'End.dat' ));

Region1LB01 = -200;
Region1UB01 = -100;
Region2UB01 = 10;
Region3UB01 = 80;
Region4UB01 = 150;
Region5UB01 = 200;

Region1LBi01 = round(((Region1LB01 - scenbeg)/dx)) + 1;
Region1UBi01 = round(((Region1UB01 - scenbeg)/dx)) + 1;
Region2UBi01  = round(((Region2UB01 - scenbeg)/dx)) + 1;
Region3UBi01 = round(((Region3UB01  - scenbeg)/dx)) + 1;
Region4UBi01 = round(((Region4UB01  - scenbeg)/dx)) + 1;
Region5UBi01 = round(((Region5UB01  - scenbeg)/dx)) + 1;

t = End(1,1);
x01 = [End(Region1LBi01:lessdenselinesep:Region1UBi01,2); ...
    End(Region1UBi01:denselinesep:Region2UBi01,2); ...
    End(Region2UBi01:lessdenselinesep:Region3UBi01,2); ...
    End(Region3UBi01:denselinesep:Region4UBi01,2); ...
    End(Region4UBi01:lessdenselinesep:Region5UBi01,2)];
    
h01 = [End(Region1LBi01:lessdenselinesep:Region1UBi01,3); ...
    End(Region1UBi01:denselinesep:Region2UBi01,3); ...
    End(Region2UBi01:lessdenselinesep:Region3UBi01,3); ...
    End(Region3UBi01:denselinesep:Region4UBi01,3); ...
    End(Region4UBi01:lessdenselinesep:Region5UBi01,3)];

G01 = [End(Region1LBi01:lessdenselinesep:Region1UBi01,4); ...
    End(Region1UBi01:denselinesep:Region2UBi01,4); ...
    End(Region2UBi01:lessdenselinesep:Region3UBi01,4); ...
    End(Region3UBi01:denselinesep:Region4UBi01,4); ...
    End(Region4UBi01:lessdenselinesep:Region5UBi01,4)];

u01 = [End(Region1LBi01:lessdenselinesep:Region1UBi01,5); ...
    End(Region1UBi01:denselinesep:Region2UBi01,5); ...
    End(Region2UBi01:lessdenselinesep:Region3UBi01,5); ...
    End(Region3UBi01:denselinesep:Region4UBi01,5); ...
    End(Region4UBi01:lessdenselinesep:Region5UBi01,5)];


wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/Examples/DB/03/";
End = importdata(strcat(wdir,'End.dat' ));

Region1LB02 = -200;
Region1UB02 = 200;

Region1LBi02 = round(((Region1LB02 - scenbeg)/dx)) + 1;
Region1UBi02 = round(((Region1UB02 - scenbeg)/dx)) + 1;

t = End(1,1);
x02 = End(Region1LBi02:lessdenselinesep:Region1UBi02,2);
h02 = End(Region1LBi02:lessdenselinesep:Region1UBi02,3);
G02 = End(Region1LBi02:lessdenselinesep:Region1UBi02,4);
u02 = End(Region1LBi02:lessdenselinesep:Region1UBi02,5);

wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/Examples/DB/02/";
End = importdata(strcat(wdir,'End.dat' ));

Region1LB03 = -200;
Region1UB03 = 200;

Region1LBi03 = round(((Region1LB03 - scenbeg)/dx)) + 1;
Region1UBi03 = round(((Region1UB03 - scenbeg)/dx)) + 1;

t = End(1,1);
x03 = End(Region1LBi03:lessdenselinesep:Region1UBi03,2);
h03 = End(Region1LBi03:lessdenselinesep:Region1UBi03,3);
G03 = End(Region1LBi03:lessdenselinesep:Region1UBi03,4);
u03 = End(Region1LBi03:lessdenselinesep:Region1UBi03,5);

figure;
plot(x00,h00,'-b');
hold on;
plot(x01,h01,'-r');
plot(x02,h02,'-g');
plot(x03,h03,'-k');
xlabel('x (m)');
ylabel('h (m)');
axis([-200 200 1.0,2.2]);
xticks([-200,-150,-100,-50,0,50,100,150,200]);
yticks([1,1.2,1.4,1.6,1.8,2,2.2]);
legend('hide')
matlab2tikz('hEx.tex');
close all;

figure;
plot(x00,G00,'-b');
hold on;
plot(x01,G01,'-r');
plot(x02,G02,'-g');
plot(x03,G03,'-k');
xlabel('x (m)');
ylabel('G (m^2/s)');
axis([-200 200 -15 15]);
xticks([-200,-150,-100,-50,0,50,100,150,200]);
yticks([-15,-10,-5,0,5,10,15]);
legend('hide')
matlab2tikz('GEx.tex');
close all;


figure;
plot(x00,u00,'-b');
hold on;
plot(x01,u01,'-r');
plot(x02,u02,'-g');
plot(x03,u03,'-k');
xlabel('x (m)');
ylabel('u (m/s)');
axis([-200 200 -0.5 2.5]);
xticks([-200,-150,-100,-50,0,50,100,150,200]);
yticks([-0.5,0,0.5,1,1.5,2,2.5]);
legend('hide')
matlab2tikz('uEx.tex');
close all;

