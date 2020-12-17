% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/AnalyticSolutions/SerreSoliton/";

a0 = 1.0;
a1 = 0.7;
g = 9.81;
c = sqrt(g*(a0 + a1));
htop = a0 + a1;
utop = c*(1 - a0/htop);
Gtop = utop*htop;

n = 12;
peakherror = zeros(n+1,1);
peakuerror = zeros(n+1,1);
peakGerror = zeros(n+1,1);
dxs = zeros(n+1,1);

for k = 0:n
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

    PeakLochj = find(h==max(h));
    PeakLocGj = find(G==max(G));
    PeakLocuj = find(u==max(u));
    hj = h(PeakLochj);
    Gj = G(PeakLocGj);
    uj = u(PeakLocuj);
    peakh = hj;
    peakG = Gj;
    peaku = uj;
    peakherror(k+1) = sqrt((peakh - htop)^2/(htop)^2);
    peakuerror(k+1) = sqrt((peaku - utop)^2/(utop)^2);
    peakGerror(k+1) = sqrt((peakG - Gtop)^2/(Gtop)^2);
    dxs(k+1) = dx;


end

loglog(dxs,peakherror,'s b', dxs,peakuerror,'^ k' );
xlabel('\Delta x')
ylabel('Peak Error')
legend('hide')
axis([10^-4 1 10^-9 10]);
xticks([10^-4,10^-3,10^-2,10^-1,10^0,10]);
yticks([10^-9,10^-7,10^-5,10^-3,10^-1,10]);
 matlab2tikz('PeakError.tex');

