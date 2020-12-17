% Process Fortran Outputs

clc;
clear all;
close all;

% Get list of directories to loop over when reading data
wdir = "/home/jp/Documents/Work/PostDoc/Projects/Steve/1DWaves/RegularisedSerre/Data/RAW/Presentations/CTAC/Animation/SerreSoliton/00/";


kstart = 1;
ksep = 10;
kend = 3001;

a0 = 1.0;
a1 = 0.7;
g = 9.81;


sep =1;

k = kstart;
SolFile = strcat(wdir,compose("%8d",k),'.dat' );
SolT = importdata(SolFile);
tcurr = SolT(1,1);
t = tcurr;
x = SolT(1:sep:end,2);
H = SolT(1:sep:end,3);
HA = hsoliton(x,tcurr,a0,a1,g);

for k = kstart + 1:ksep:kend
    SolFile = strcat(wdir,compose("%8d",k),'.dat' );

    SolT = importdata(SolFile);
    tcurr = SolT(1,1);
    t = [t,tcurr];
    H = [H,SolT(1:sep:end,3)];
    HA = [HA,hsoliton(x,tcurr,a0,a1,g)];
    
end

fontsizenum = 14;
sizeH = size(H);
sizeHA = size(HA);
Vid = VideoWriter('Soliton');
open(Vid);
figure('units','pixels','position',[0 0 1280 720]) ;
for i = 1: sizeH(2)
    HAcurr = HA(:,i);
    Hcurr = H(:,i);
    area(x,HAcurr,'FaceColor', 'b', 'FaceAlpha',0.3,'EdgeColor','b','LineStyle','-' , 'LineWidth',1);
    axis([-50,150,0.8,1.8]);
    xticks([-50,0,50,100,150]);
    yticks([0.8,1,1.2,1.4,1.6,1.8]);
    xlabel('x (m)')
    ylabel('h (m)')
    text(-40,1.65,strcat('t = ',sprintf("%5.2f",t(i))),'fontsize',fontsizenum,'fontname','cmr12');
    hold on;
    plot(x,Hcurr,':r' ,'LineWidth',5);
    set(gca,'fontsize', fontsizenum,'fontname','cmr12');
    hold off;
    frame = getframe(gcf);
    writeVideo(Vid,frame);
    clf ;

end
close(Vid);


function ha = hsoliton(x,t,a0,a1,g)
    c =  sqrt(g*(a0 + a1));
    k =  sqrt(3*a1) / (2*a0*sqrt(a0 + a1));
    Phi = x - c*t;
    ha = a0 + a1*sech(k*Phi).^2;
end


