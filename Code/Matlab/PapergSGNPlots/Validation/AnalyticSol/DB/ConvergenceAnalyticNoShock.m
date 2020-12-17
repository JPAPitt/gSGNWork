% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/AnalyticSolutions/DBSWWETheta1/";

h0 = 2.0;
h1 = 1.0;
g = 9.81;
nBC = 6;

n = 12;
L2hDown = zeros(n+1,1);
L2uDown = zeros(n+1,1);
L2GDown = zeros(n+1,1);

L2hMiddle = zeros(n+1,1);
L2uMiddle = zeros(n+1,1);
L2GMiddle = zeros(n+1,1);

dxs = zeros(n+1,1);

func1 = @(x) h2DB(x,h0,h1,g);
h2 = fzero(func1,h0);
u2 = 2*(sqrt(g*h0) -sqrt(g*h2) );
S2 = 2*h2/(h2 - h1)*(sqrt(g*h0) - sqrt(g*h2));

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
    
    hA = zeros(size(h));
    uA = zeros(size(h));
    
    for i = 1:length(h)
        if (x(i) <= -t*sqrt(g*h0))
           hA(i) = h0;
           uA(i) = 0;
        elseif ((x(i) > -t*sqrt(g*h0)) & (x(i) < t*(u2 - sqrt(g*h2))))
           hA(i) = 4.0/(9*g) *(sqrt(g*h0) - x(i)/(2*t))^2;
           uph = 2.0/3.0*(sqrt(g*h0)*(x(i) + 0.5*dx) + 1.0/2.0*(x(i) + 0.5*dx)^2/t);
           umh = 2.0/3.0*(sqrt(g*h0)*(x(i) - 0.5*dx) + 1.0/2.0*(x(i) - 0.5*dx)^2/t);
           uA(i) = (uph -umh) / dx;
        elseif ((x(i) >= t*(u2 - sqrt(g*h2)) ) & (x(i) <= t*S2))
           hA(i) = h2  ;  
           uA(i) = u2 ; 
        else
           hA(i) = h1;
           uA(i) = 0;
        end
    end
    
    GA = uA.*hA;
    
    Tol = 30;
    %Middle Region
    MiddleLB = t*(u2 - sqrt(g*h2)) + Tol;
    MiddleUB = t*S2 - Tol;
    
    MiddleLBindex = find( x > MiddleLB, 1 );
    
    MiddleUBindex = find( x > MiddleUB, 1 );
    
    hMid = h(MiddleLBindex: MiddleUBindex);
    hAMid = hA(MiddleLBindex: MiddleUBindex);
    
    uMid = u(MiddleLBindex: MiddleUBindex);
    uAMid = uA(MiddleLBindex: MiddleUBindex);
    
    GMid = G(MiddleLBindex: MiddleUBindex);
    GAMid = GA(MiddleLBindex: MiddleUBindex);
   
%     
    L2hMiddle(k+1) = norm(hAMid - hMid,2)/ norm(hAMid,2);
    L2uMiddle(k+1) = norm(uAMid - uMid,2)/ norm(uAMid,2);
    L2GMiddle(k+1) = norm(GAMid - GMid,2)/ norm(GAMid,2);
    dxs(k+1) = dx;
    
    %Down Region
    DownLB = -t*sqrt(g*h0) + Tol;
    DownUB = t*(u2 - sqrt(g*h2))- Tol;
    
    DownLBindex = find( x > DownLB, 1 );
    
    DownUBindex = find( x > DownUB, 1 );
    
    hDown = h(DownLBindex: DownUBindex);
    hADown = hA(DownLBindex: DownUBindex);
    
    uDown = u(DownLBindex: DownUBindex);
    uADown = uA(DownLBindex: DownUBindex);
    
    GDown = G(DownLBindex: DownUBindex);
    GADown = GA(DownLBindex: DownUBindex);
   
%     
    L2hDown(k+1) = norm(hADown - hDown,2)/ norm(hADown,2);
    L2uDown(k+1) = norm(uADown - uDown,2)/ norm(uADown,2);
    L2GDown(k+1) = norm(GADown - GDown,2)/ norm(GADown,2);
    dxs(k+1) = dx;



end
% 
% figure;
% loglog(dxs,L2hMiddle,'s b',dxs,L2uMiddle,'^ r', dxs,L2GMiddle,'o k');
% grid off;
% legend('hide');
% xlabel('\Delta x');
% ylabel('L_2');
% axis([10^-3 1 10^-10 10^-4]);
% xticks([10^-4,10^-3,10^-2,10^-1,10^0,10]);
% yticks([10^-10,10^-09,10^-08,10^-07, 10^-06,10^-05,10^-4]);
% matlab2tikz('ConvergenceInConstantState.tex');



function y = h2DB(x,h0,h1,g)
y = x - h1/2*( sqrt( 1+ 8*(2*x/(x - h1) *((sqrt(g*h0) - sqrt(g*x))/ sqrt(g*h1) ))^2 ) -1) ;
end

