% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdirbase = "H:\Work\UoA\Github\gSGNWork\Results";
ExpRes = "\Data\Segur-1cm\";
Sim1Res = "\Experimental\SegurSGN\1cm\01\";
Sim2Res = "\Experimental\SeguriSGN\1cm\01\";
Sim3Res = "\Experimental\SegurSWWE\1cm\01\";

path(pathdef)
path(path,'..')
path(path,'H:\Work\UoA\Github\SingleFloeLinearModel\CommonFunctions\Mat2Tikzsrc')


x0s = [0,5,10,15,20];
letters = ["a","b","c","d","e"];


for i = 1:5
    ExpFig = importdata(strcat(wdirbase,ExpRes,'fig2',letters(i),'.txt' ));
    Sim1Loc = importdata(strcat(wdirbase,Sim1Res,"TS_ ",num2str(i),'.dat' ));
    Sim2Loc = importdata(strcat(wdirbase,Sim2Res,"TS_ ",num2str(i),'.dat' ));
    Sim3Loc = importdata(strcat(wdirbase,Sim3Res,"TS_ ",num2str(i),'.dat' ));



    x0 = x0s(i);
    g = 9.81;
    h0 = 0.1;

    Time1 = Sim1Loc(:,1);
    h = Sim1Loc(:,3);

    ModTime1 = Time1*sqrt(g/h0) - x0/h0;
    Modh1 =  (h - h0) / h0;

    Time1 = Sim2Loc(:,1);
    h = Sim2Loc(:,3);

    ModTime2 = Time1*sqrt(g/h0) - x0/h0;
    Modh2 =  (h - h0) / h0;


    Time1 = Sim3Loc(:,1);
    h = Sim3Loc(:,3);
    
    ModTime3 = Time1*sqrt(g/h0) - x0/h0;
    Modh3 = (h - h0) / h0;

    figure();
    hold on;
    plot(ModTime1,Modh1,'-b');
    plot(ModTime2,Modh2,'-k');
    plot(ModTime3,Modh3,'--g');
    plot(ExpFig(:,1),ExpFig(:,2) / (1.5),'-r');
    xlabel('t')
    ylabel('h')
    xlim([0,200])
    ylim([-0.25,0.15])

    str1 = ['Seg1cm',char(letters(i)),'.tex'];
    cleanfigure;
    matlab2tikz(str1); 
end

