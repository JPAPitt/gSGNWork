% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Presentations/CTAC/Animation/DB6/00/";


kstart = 1;
ksep = 2;
kend = 2001;


sep =1;

k = kstart;
SolFile = strcat(wdir,compose("%8d",k),'.dat' );
SolT = importdata(SolFile);
tcurr = SolT(1,1);
t = tcurr;
x = SolT(1:sep:end,4);
H = SolT(1:sep:end,5);
beta1 = SolT(1,2) + 2.0/3.0;
beta2 = SolT(1,3);

for k = kstart + 1:ksep:kend
    SolFile = strcat(wdir,compose("%8d",k),'.dat' );

    SolT = importdata(SolFile);
    tcurr = SolT(1,1);
    beta1curr =SolT(1,2) + 2.0/3.0;
    beta2curr =SolT(1,3) ;
    beta1 = [beta1,beta1curr];
    beta2 = [beta2,beta2curr];
    t = [t,tcurr];
    H = [H,SolT(1:sep:end,5)];
    
end

fontsizenum = 14;
sizeH = size(H);
Vid = VideoWriter('Dambreak');
open(Vid);
figure('units','pixels','position',[0 0 1280 720]) ;
for i = 1: sizeH(2)
    Hcurr = H(:,i);
    area(x,Hcurr,'FaceColor', 'b', 'FaceAlpha',0.3,'EdgeColor','b','LineStyle','-' , 'LineWidth',1);
    axis([-1000,1000,0.8,2.2]);
    xticks([-1000,-500,0,500,1000]);
    yticks([0.8,1,1.2,1.4,1.6,1.8,2,2.2]);
    xlabel('x (m)')
    ylabel('h (m)')
    text(-900,2.1,strcat('t = ',sprintf("%5.2f",t(i))),'fontsize',fontsizenum,'fontname','cmr12');
    text(-500,2.1,strcat('\beta_1 = ',sprintf("%5.2f",beta1(i))),'fontsize',fontsizenum,'fontname','cmr12');
    text(-100,2.1,strcat('\beta_2 = ',sprintf("%5.2f",beta2(i))),'fontsize',fontsizenum,'fontname','cmr12');
    hold on;
    set(gca,'fontsize', fontsizenum,'fontname','cmr12');
    hold off;
    frame = getframe(gcf);
    writeVideo(Vid,frame);
    clf ;
end
close(Vid);




