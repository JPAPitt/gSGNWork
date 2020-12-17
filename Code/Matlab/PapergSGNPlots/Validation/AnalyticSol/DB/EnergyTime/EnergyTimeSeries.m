% Process Fortran Outputs

clc;
clear all;
close all;


kstart = 0;
ksep = 2;
kend = 12;
figure;
hold on;
for k = kstart:ksep:kend
    NumStr = compose("%2.2d",k);
    TimeSeries = importdata(strcat('DB_DX',NumStr,'.dat'));
    t = TimeSeries(:,1);
    H = TimeSeries(:,6);
    dx= TimeSeries(1,2);
    plot(t,H,'.','DisplayName',strcat('dx = ',num2str(dx)));
    
end

DBEs = zeros(1,length(t));
for tv = 1:length(t)
    DBEs(tv) = DBEN(t(tv),2,1,9.81);
end

plot(t,DBEs,'--k','DisplayName',strcat('Analytic Solution'))

xlabel('t (s)');
ylabel('H (m^3/s^2)');
legend('show');
axis([0 35 6105 6140]);
xticks([0,5,10,15,20,25,30,35]);
yticks([6105,6110,6115,6120,6125,6130,6135,6140]);
legend('show');
hold off
matlab2tikz('EnergyOverTime.tex');



function h2func = f2func(x,h0,h1,g)
    h2func = x - h1/2.0*(sqrt(1 + 8*( ((2*x)/(x-h1))*(( sqrt(g*h0) - sqrt(g*x))/ (sqrt(g*h1))) )^2) -1); 
end

function DBE = DBEN(t,h0,h1,g)

    if t == 0
        ghEp1 = h0^2*(0 - -250);
        ghEp2 = h1^2*(250);
        DBE = g*(ghEp1 + ghEp2)/2;
    else
        fun = @(x) f2func(x,h0,h1,g);
        h2 = fzero(fun,h0) ;   
        u2 = 2*(sqrt(g*h0) - sqrt(g*h2));
        S2 = (2*h2/ (h2 - h1))*(sqrt(g*h0) - sqrt(g*h2));

        x0 = -250;
        x1 = -t*sqrt(g*h0);
        x2 = t*(u2 - sqrt(g*h2));
        x3 = t*S2;
        x4 = 250;

        %gh^2
        ghEp1 = h0^2*(x1 - x0);
        Int2 = -2.0/5.0*t*(sqrt(g*h0) - x2/(2*t))^5;
        Int1 = -2.0/5.0*t*(sqrt(g*h0) - x1/(2*t))^5;
        ghEp2 = (4/(9*g))^2*(Int2 - Int1);
        ghEp3 = h2^2*(x3 - x2);
        ghEp4 = h1^2*(x4 - x3);

        gh = g*(ghEp1 + ghEp2 + ghEp3 + ghEp4);

        %uh^2
        Int2 = (4*g^2*h0^2*t^4*x2 - g*h0*t^2*x2^2*(x2 - 2*t*sqrt(g*h0)) +x2^4/10.0*(2*x2 - 5*t*sqrt(g*h0))) / (4*t^4);
        Int1 = (4*g^2*h0^2*t^4*x1 - g*h0*t^2*x1^2*(x1 - 2*t*sqrt(g*h0)) +x1^4/10.0*(2*x1 - 5*t*sqrt(g*h0))) / (4*t^4);
        u2hEp2 =  (2.0/3.0)^2*(4/(9*g))*(Int2 - Int1);
        uhEp3 = u2^2*h2*(x3 - x2);

        uh = u2hEp2+uhEp3;
        DBE = (gh + uh)/2;
    end
end

    
