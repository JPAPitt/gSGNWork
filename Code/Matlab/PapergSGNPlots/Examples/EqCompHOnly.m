% Process Fortran Outputs

clc;
clear all;
close all;

wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/Examples/DBNew/00/";

param = fileread(strcat(wdir ,'Params.dat' ));
dxstr = extractBetween(param,"dx","tstart");
dx = str2double(dxstr{1,1});

scenbeg = -250;
scenend = 250;


Region1LB00 = -250;
Region1UB00 = -100;
Region2UB00 = 150;
Region3UB00 = 250;


denselinesep = 1;
lessdenselinesep = 5;

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

wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/Examples/DBNew/01/";
End = importdata(strcat(wdir,'End.dat' ));

Region1LB01 = -250;
Region1UB01 = -100;
Region2UB01 = 10;
Region3UB01 = 80;
Region4UB01 = 150;
Region5UB01 = 250;

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


wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/Examples/DBNew/03/";
End = importdata(strcat(wdir,'End.dat' ));

Region1LB02 = -250;
Region1UB02 = 250;

Region1LBi02 = round(((Region1LB02 - scenbeg)/dx)) + 1;
Region1UBi02 = round(((Region1UB02 - scenbeg)/dx)) + 1;

t = End(1,1);
x02 = End(Region1LBi02:lessdenselinesep:Region1UBi02,2);
h02 = End(Region1LBi02:lessdenselinesep:Region1UBi02,3);
G02 = End(Region1LBi02:lessdenselinesep:Region1UBi02,4);
u02 = End(Region1LBi02:lessdenselinesep:Region1UBi02,5);

wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/Examples/DBNew/02/";
End = importdata(strcat(wdir,'End.dat' ));

Region1LB03 = -250;
Region1UB03 = 250;

Region1LBi03 = round(((Region1LB03 - scenbeg)/dx)) + 1;
Region1UBi03 = round(((Region1UB03 - scenbeg)/dx)) + 1;

t = End(1,1);
x03 = End(Region1LBi03:lessdenselinesep:Region1UBi03,2);
h03 = End(Region1LBi03:lessdenselinesep:Region1UBi03,3);
G03 = End(Region1LBi03:lessdenselinesep:Region1UBi03,4);
u03 = End(Region1LBi03:lessdenselinesep:Region1UBi03,5);

wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/Examples/DBNew/04/";
End = importdata(strcat(wdir,'End.dat' ));

Region1LB04 = -250;
Region1UB04 = -234;
Region2UB04 = -160;
Region3UB04 = 140;
Region4UB04 = 170;
Region5UB04 = 250;

Region1LBi04 = round(((Region1LB04 - scenbeg)/dx)) + 1;
Region1UBi04 = round(((Region1UB04 - scenbeg)/dx)) + 1;
Region2UBi04  = round(((Region2UB04 - scenbeg)/dx)) + 1;
Region3UBi04 = round(((Region3UB04  - scenbeg)/dx)) + 1;
Region4UBi04 = round(((Region4UB04  - scenbeg)/dx)) + 1;
Region5UBi04 = round(((Region5UB04  - scenbeg)/dx)) + 1;

t = End(1,1);
x04 = [End(Region1LBi04:lessdenselinesep:Region1UBi04,2); ...
    End(Region1UBi04:denselinesep:Region2UBi04,2); ...
    End(Region2UBi04:lessdenselinesep:Region3UBi04,2); ...
    End(Region3UBi04:denselinesep:Region4UBi04,2); ...
    End(Region4UBi04:lessdenselinesep:Region5UBi04,2)];
    
h04 = [End(Region1LBi04:lessdenselinesep:Region1UBi04,3); ...
    End(Region1UBi04:denselinesep:Region2UBi04,3); ...
    End(Region2UBi04:lessdenselinesep:Region3UBi04,3); ...
    End(Region3UBi04:denselinesep:Region4UBi04,3); ...
    End(Region4UBi04:lessdenselinesep:Region5UBi04,3)];

G04 = [End(Region1LBi04:lessdenselinesep:Region1UBi04,4); ...
    End(Region1UBi04:denselinesep:Region2UBi04,4); ...
    End(Region2UBi04:lessdenselinesep:Region3UBi04,4); ...
    End(Region3UBi04:denselinesep:Region4UBi04,4); ...
    End(Region4UBi04:lessdenselinesep:Region5UBi04,4)];

u04 = [End(Region1LBi04:lessdenselinesep:Region1UBi04,5); ...
    End(Region1UBi04:denselinesep:Region2UBi04,5); ...
    End(Region2UBi04:lessdenselinesep:Region3UBi04,5); ...
    End(Region3UBi04:denselinesep:Region4UBi04,5); ...
    End(Region4UBi04:lessdenselinesep:Region5UBi04,5)];

figure;
plot(x04,h04,'-b');
xlabel('x (m)');
ylabel('h (m)');
axis([-250 250 0.8,2.2]);
xticks([-250,-125,0,125,250]);
yticks([0.8,1,1.2,1.4,1.6,1.8,2,2.2]);
legend('hide')
matlab2tikz('hEx04.tex');

