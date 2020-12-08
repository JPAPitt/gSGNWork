% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
%wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGNForcedLimAll/ConstantBeta/AnaSolSolitonLoop/06/";
 wdir = "./FD2C/Validation/AnalyticSolutions/SerreSoliton/06/"

End = importdata(strcat(wdir,'End.dat' ));
EndA = importdata(strcat(wdir,'EndAna.dat' ));


t = End(1,1);
x = End(:,2);
h = End(:,3);
u = End(:,4);

t = EndA(1,1);
xA = EndA(:,2);
hA = EndA(:,3);
uA = EndA(:,4);




figure;
plot(xA,hA,'-b');
hold on;
plot(x,h,'.r');
hold off

figure;
plot(xA,uA,'-b');
hold on;
plot(x,u,'.r');
hold off









