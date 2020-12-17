% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/AnalyticSolutions/DBSWWETheta1/";

h0 = 2.0;
h1 = 1.0;
g = 9.81;

LocationUB = 100;
LocationLB = 0;

n = 12;
Stateherror = zeros(n+1,1);
Stateuerror = zeros(n+1,1);
dxs = zeros(n+1,1);

func1 = @(x) h2DB(x,h0,h1,g);
h2 = fzero(func1,h0);
u2 = 2*(sqrt(g*h0) -sqrt(g*h2) );
% h2 = fzero(func1,h1);

for k = 0:12
    dxNumStr = compose("%2.2d",k);
    expdir = strcat(wdir,dxNumStr,'/');
    EndF = importdata(strcat(expdir, 'End.dat'));
    t = EndF(1,1);
    x = EndF(:,2);
    h = EndF(:,3);
    G = EndF(:,4);
    u = EndF(:,5);
       
    param = fileread(strcat(expdir ,'Params.dat' ));
    dxstr = extractBetween(param,"dx","tstart");
    dx = str2double(dxstr{1,1});
    
    xLB = find( x > LocationLB -0.5*dx, 1 );
    xUB = find( x > LocationUB -0.5*dx , 1 );
    
    meanh = mean(h(xLB:xUB));
    meanu = mean(u(xLB:xUB));
%     h2 = fzero(func1,meanh);
    Stateherror(k+1) = abs(meanh - h2)/(h2);
    Stateuerror(k+1) = abs(meanu - u2)/(u2);
    dxs(k+1) = dx;



end

loglog(dxs,Stateherror,'s b',dxs,Stateuerror,'^ k', dxs, 0.0001*dxs,'-', dxs, 0.0001*dxs.^2,'-');
legend('h', 'u', '1st', '2nd');

function y = h2DB(x,h0,h1,g)
y = x - h1/2*( sqrt( 1+ 8*(2*x/(x - h1) *((sqrt(g*h0) - sqrt(g*x))/ sqrt(g*h1) ))^2 ) -1) ;
end


