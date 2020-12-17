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
SpeedUB = zeros(n+1,1);
SpeedLB = zeros(n+1,1);
dxs = zeros(n+1,1);

func1 = @(x) h2DB(x,h0,h1,g);
h2 = fzero(func1,h0);
u2 = 2*(sqrt(g*h0) -sqrt(g*h2) );
S2 = 2*h2/(h2 - h1)*(sqrt(g*h0) - sqrt(g*h2));
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
    
    UBloc = find( h < h1  + (h2 - h1)/10.0, 1 );
    LBloc = find( h < h1  + 9*(h2 - h1)/10.0, 1 );
    xubloc = x(UBloc);
    xlbloc = x(LBloc);
    
    SpeedUB(k+1) = (xubloc / t)/S2;
    SpeedLB(k+1) = (xlbloc / t)/S2;
    
    dxs(k+1) = dx;



end

semilogx(dxs,SpeedUB,'x b',dxs,SpeedLB,'x r');
hold on;
plot([10^(-4),10],[1,1],'--k');
legend('hide');
xlabel('\Delta x');
ylabel('S_2 / S_2 ');
axis([10^-3 10 0.96 1.08]);
xticks([10^-3,10^-2,10^-1,10^0,10]);
yticks([0.96, 0.98,1,1.02,1.04,1.06,1.08]);
matlab2tikz('ShockSpeedEstimates.tex');
close all;


function y = h2DB(x,h0,h1,g)
y = x - h1/2*( sqrt( 1+ 8*(2*x/(x - h1) *((sqrt(g*h0) - sqrt(g*x))/ sqrt(g*h1) ))^2 ) -1) ;
end


