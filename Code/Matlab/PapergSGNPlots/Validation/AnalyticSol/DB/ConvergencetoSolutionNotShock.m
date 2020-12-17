% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Models/gSGN_NoLim/AnalyticSolutions/DBSWWE/";

h0 = 2.0;
h1 = 1.0;
g = 9.81;
nBC = 6;

n = 12;
L2h = zeros(n+1,1);
dxs = zeros(n+1,1);

HR_k = 12;
dxNumStr = compose("%2.2d",HR_k);
expdir = strcat(wdir,dxNumStr,'/');
EndF = importdata(strcat(expdir, 'End.dat'));
t = EndF(1,1);
xHR = EndF(:,2);
hHR = EndF(:,3);
GHR = EndF(:,4);
uHR = EndF(:,5);


for k = 0:1
    
    dxNumStr = compose("%2.2d",(HR_k -1) - k);
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
    
    newlen = length(xHR) - 2*nBC;
    hlowonhigh = zeros(newlen,1);
    
    Current_LowxCell = nBC + 1;
    for j = 1: length(hlowonhigh)
        
        if xHR(j + nBC) > x(Current_LowxCell) + 0.5*dx
            Current_LowxCell = Current_LowxCell + 1;
        end        
        hlowonhigh(j) = h(Current_LowxCell) + (xHR(j + nBC) -x(Current_LowxCell) )*(h(Current_LowxCell + 1) - h(Current_LowxCell - 1))/(2*dx);
        
    end
    
    nx = xHR(nBC +1 : end- nBC);
    nh = hHR(nBC +1 : end- nBC);
%     
    L2h(k+1) = norm(hlowonhigh - nh,2)/ norm(nh,2);
    dxs(k+1) = dx;



end
% 
loglog(dxs,L2h,'s b', dxs, 0.0001*dxs,'-', dxs, 0.0001*dxs.^2,'-');
legend('h', '1st', '2nd');



