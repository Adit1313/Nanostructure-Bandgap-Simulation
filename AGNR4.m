close all;
clear all;
clc;

n = 4;
N_Ring_Unit = 2;

N_Total_Plot = (n*N_Ring_Unit*3);

a0 = 1.42;
d1=a0;
d2=2.45951214;
d3=2.84;
mrg = 0.01;

eps=0;
t1=-2.65;
t2=0;
t3=-0.27;

D_w = (n-1)*(a0*sqrt(3)/2);
Unit_L = N_Ring_Unit*1.5*a0;
fprintf('Ribbon Width = %.2f nm\n\n',D_w/10);

atom_x=zeros(1,N_Total_Plot);
atom_y=zeros(1,N_Total_Plot);

atom_x(1)=0;     atom_x(2)=d1/2;           atom_x(3)=0;            atom_x(4)=d1/2;
atom_y(1)=0;     atom_y(2)=d1*sqrt(3)/2;   atom_y(3)=d1*sqrt(3);   atom_y(4)=1.5*d1*sqrt(3);

atom_x(5)=2*d1;  atom_x(6)=d1*1.5;         atom_x(7)=2*d1;         atom_x(8)=d1*1.5;
atom_y(5)=0;     atom_y(6)=d1*sqrt(3)/2;   atom_y(7)=d1*sqrt(3);   atom_y(8)=1.5*d1*sqrt(3);

atom_x(9)=0+3*d1;      atom_x(10)=d1/2+3*d1;     atom_x(11)=0+3*d1;      atom_x(12)=d1/2+3*d1;
atom_y(9)=0;           atom_y(10)=d1*sqrt(3)/2;  atom_y(11)=d1*sqrt(3);  atom_y(12)=1.5*d1*sqrt(3);

atom_x(13)=2*d1+3*d1;  atom_x(14)=d1*1.5+3*d1;   atom_x(15)=2*d1+3*d1;   atom_x(16)=d1*1.5+3*d1;
atom_y(13)=0;          atom_y(14)=d1*sqrt(3)/2;  atom_y(15)=d1*sqrt(3);  atom_y(16)=1.5*d1*sqrt(3);

atom_x(17)=0+6*d1;     atom_x(18)=d1/2+6*d1;     atom_x(19)=0+6*d1;      atom_x(20)=d1/2+6*d1;
atom_y(17)=0;          atom_y(18)=d1*sqrt(3)/2;  atom_y(19)=d1*sqrt(3);  atom_y(20)=1.5*d1*sqrt(3);

atom_x(21)=2*d1+6*d1;  atom_x(22)=d1*1.5+6*d1;   atom_x(23)=2*d1+6*d1;   atom_x(24)=d1*1.5+6*d1;
atom_y(21)=0;          atom_y(22)=d1*sqrt(3)/2;  atom_y(23)=d1*sqrt(3);  atom_y(24)=1.5*d1*sqrt(3);

for i=1:N_Total_Plot
    plot(atom_x(i),atom_y(i),'o','MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor',[0,0,1]);hold on;
    axis equal;
end
title('4-AGNR');
line([(Unit_L - 0.5*d1) (Unit_L - 0.5*d1)],[-1.5 (D_w+1.5)],'Color','k','LineWidth',1);
line([(2*Unit_L - 0.5*d1) (2*Unit_L - 0.5*d1)],[-1.5 (D_w+1.5)],'Color','k','LineWidth',1);

fprintf('Making Hamiltonians...');

N_Unit = n*N_Ring_Unit;

atom_x0 = zeros(1,N_Unit);
atom_y0 = zeros(1,N_Unit);

atom_xl = zeros(1,N_Unit);
atom_yl = zeros(1,N_Unit);

atom_xr = zeros(1,N_Unit);
atom_yr = zeros(1,N_Unit);

H0=zeros(N_Unit,N_Unit);
HL=zeros(N_Unit,N_Unit);
HR=zeros(N_Unit,N_Unit);

for i=1:N_Unit
    atom_xl(i)=atom_x(i);
    atom_yl(i)=atom_y(i);
end

for i=1:N_Unit
    atom_x0(i)=atom_x(i+N_Unit);
    atom_y0(i)=atom_y(i+N_Unit);
end

for i=1:N_Unit
    atom_xr(i)=atom_x(i+2*N_Unit);
    atom_yr(i)=atom_y(i+2*N_Unit);
end

for j=1:N_Unit
    for i=1:N_Unit
        
        dis = sqrt( (atom_x0(j)-atom_x0(i))^2 + (atom_y0(j)-atom_y0(i))^2 );
        
        if  dis < mrg
            H0(j,i) = eps;
            
        elseif (d1-mrg<dis) && (dis<d1+mrg)
            H0(j,i) = t1;
            
        elseif (d2-mrg<dis) && (dis<d2+mrg)
            H0(j,i) = t2;
            
        elseif (d3-mrg<dis) && (dis<d3+mrg)
            H0(j,i) = t3;
            
        else
            H0(j,i) = 0;
        end
    end
end

for j=1:N_Unit
    for i=1:N_Unit
        
        dis = sqrt( (atom_x0(j)-atom_xl(i))^2 + (atom_y0(j)-atom_yl(i))^2 );
        
        if  dis < mrg
            HL(j,i) = eps;
            
        elseif (d1-mrg<dis) && (dis<d1+mrg)
            HL(j,i) = t1;
            
        elseif (d2-mrg<dis) && (dis<d2+mrg)
            HL(j,i) = t2;
            
        elseif (d3-mrg<dis) && (dis<d3+mrg)
            HL(j,i) = t3;
            
        else
            HL(j,i) = 0;
        end
    end
end

for j=1:N_Unit
    for i=1:N_Unit
        
        dis = sqrt( (atom_x0(j)-atom_xr(i))^2 + (atom_y0(j)-atom_yr(i))^2 );
        
        if  dis < mrg
            HR(j,i) = eps;
            
        elseif (d1-mrg<dis) && (dis<d1+mrg)
            HR(j,i) = t1;
            
        elseif (d2-mrg<dis) && (dis<d2+mrg)
            HR(j,i) = t2;
            
        elseif (d3-mrg<dis) && (dis<d3+mrg)
            HR(j,i) = t3;
            
        else
            HR(j,i) = 0;
        end
    end
end

disp=[];
p=1;
for ka=-pi:2*pi/500:pi
    HK=H0+HR*exp(1i*ka)+HL*exp(-1i*ka);
    ee=eig(HK);
    disp=[disp,ee];
    
    for j=1:length(ee)
        if ee(j)>eps
            g(p)=ee(j);
            p=p+1;
        end
    end
end

GAP=2*(min(g)-eps); if GAP< 1e-2 GAP=0; end

figure();
plot(-pi:2*pi/500:pi,disp);
title('Band Structure of 4-AGNR');
xlabel('Ka');
ylabel('Energy(ev)');
string = sprintf('BG = %.2f eV',GAP);
legend(string);