close all;
clear;
clc;
fprintf('\t\t\tZigzag Carbon Nanotube (E-K Diagram)\n\n\n');

n = input('Input 2 * index of your ZCNT: ');  % Number of atoms on the edge of ZCNT.
N_Ring_Unit = input('Enter number of rings in ZCNT: ');     % Number of rings in unit cell

N_Total_Plot = (n*N_Ring_Unit*3);    % Number of atoms in the Plot

% Basic Distances in Carbon Honey-Comb Structure
a0 = 1.42;      % Lattice Constant (in Agnestrom)
d1=a0;          % First-neighbours distance = a0
d2=a0*sqrt(3);  % Second-neighbours distance = a0*sqrt(3)
d3=2*a0;        % Third-neighbours distance
mrg = 0.01;    % Margin needed in calculations (small value)

% Hopping Parameters
eps=0;       % On site potential for Carbon
t1=-2.65;       % First-neighbour interaction potential
t2=0;        % Second-neighbour interaction potential
t3=0;           % Third-neighbour interaction potential

D_w = ((n-1)*(a0*sqrt(3)/2));
Unit_L = N_Ring_Unit*1.5*a0;
D_r = (D_w + a0*sind(60))/(2*pi);      % Tube Radius of CNT
fprintf('Device Radius = %.2f nm \n\n',D_r/10);
fprintf('Generating atoms sites...');

% Assigning Geometric position to each atom
atom_x=zeros(1,N_Total_Plot);     % X site of atoms
atom_y=zeros(1,N_Total_Plot);     % Y site of atoms


if mod(n,2) == 0
    for j=1:(N_Ring_Unit/2)*3
        for i=1:n/2
            atom_x(2*(i-1)+1 +(j-1)*2*n)=0 + (j-1)*d1*3;       atom_x(2*i +(j-1)*2*n)=d1/2 + (j-1)*d1*3;
            atom_y(2*(i-1)+1 +(j-1)*2*n)=(i-1)*d1*sqrt(3);     atom_y(2*i +(j-1)*2*n)=(i-.5)*d1*sqrt(3);
            
            atom_x(2*(i-1)+1 +n +(j-1)*2*n)=d1*2 + (j-1)*d1*3;    atom_x(2*i +n +(j-1)*2*n)=d1*1.5 + (j-1)*d1*3;
            atom_y(2*(i-1)+1 +n +(j-1)*2*n)=(i-1)*d1*sqrt(3);     atom_y(2*i +n +(j-1)*2*n)=(i-.5)*d1*sqrt(3);
        end
    end
else
    for j=1:(N_Ring_Unit/2)*3
        for i=1:n/2
            atom_x(2*(i-1)+1 +(j-1)*2*n)=0 + (j-1)*d1*3;       atom_x(2*i +(j-1)*2*n)=d1/2 + (j-1)*d1*3;
            atom_y(2*(i-1)+1 +(j-1)*2*n)=(i-1)*d1*sqrt(3);     atom_y(2*i +(j-1)*2*n)=(i-.5)*d1*sqrt(3);
            
            atom_x(2*(i-1)+1 +n +(j-1)*2*n)=d1*2 + (j-1)*d1*3;    atom_x(2*i +n +(j-1)*2*n)=d1*1.5 + (j-1)*d1*3;
            atom_y(2*(i-1)+1 +n +(j-1)*2*n)=(i-1)*d1*sqrt(3);     atom_y(2*i +n +(j-1)*2*n)=(i-.5)*d1*sqrt(3);
        end
        
        atom_x(n +(j-1)*2*n)=0 + (j-1)*d1*3;
        atom_y(n +(j-1)*2*n)=(floor(n/2))*d1*sqrt(3);
        
        atom_x(2*n +(j-1)*2*n)=d1*2 + (j-1)*d1*3;
        atom_y(2*n +(j-1)*2*n)=(floor(n/2))*d1*sqrt(3);
    end
end

% Plotting the atoms
for i=1:N_Total_Plot
    plot(atom_x(i),atom_y(i),'o','MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor',[0,0,1]);hold on;
    axis equal;
end
title('ZCNT');

%Plot borders of unit cells
line([(Unit_L - 0.5*d1) (Unit_L - 0.5*d1)],[-1.5 (D_w+1.5)],'Color','k','LineWidth',1);
line([(2*Unit_L - 0.5*d1) (2*Unit_L - 0.5*d1)],[-1.5 (D_w+1.5)],'Color','k','LineWidth',1);

fprintf('\t\t\t-> Done!\n\n');
fprintf('Making Hamiltonians...');
N_Unit = n*N_Ring_Unit; % The number of atoms in each unit cell

atom_x0 = zeros(1,N_Unit);
atom_y0 = zeros(1,N_Unit);

atom_xl = zeros(1,N_Unit);
atom_yl = zeros(1,N_Unit);

atom_xr = zeros(1,N_Unit);
atom_yr = zeros(1,N_Unit);

% Hamiltonian matrices
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

Rows_Interval = d1*sind(60);
Roll_Over = D_w + Rows_Interval;

% Generate H0
for j=1:N_Unit
    for i=1:N_Unit
        
        dis = sqrt( (atom_x0(j)-atom_x0(i))^2 + (atom_y0(j)-atom_y0(i))^2 );
        
        %We compare the distances and assign a hopping parameter
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
    
    % Roll over atoms at the edges of AGNR to build Hamiltonian of ZCNT
    if ( (atom_y0(j) == 0) || (atom_y0(j) == Rows_Interval) )
        
        for i=1:N_Unit
            
            dis = sqrt( (atom_x0(j)-atom_x0(i))^2 + ((atom_y0(j)+Roll_Over)-atom_y0(i))^2 );
            
            if (d1-mrg<dis) && (dis<d1+mrg)
                H0(j,i) = t1;
                H0(i,j) = t1;
                
            elseif (d2-mrg<dis) && (dis<d2+mrg)
                H0(j,i) = t2;
                H0(i,j) = t2;
                
            elseif (d3-mrg<dis) && (dis<d3+mrg)
                H0(j,i) = t3;
                H0(i,j) = t3;
            end
        end
    end
end



% Generate HL
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
    
    % Roll over atoms at the edges of AGNR to biuld Hamiltonian of ZCNT
    if ( (atom_y0(j) == 0) || (atom_y0(j) == Rows_Interval) )
        
        for i=1:N_Unit
            
            dis = sqrt( (atom_x0(j)-atom_xl(i))^2 + ((atom_y0(j)+Roll_Over)-atom_yl(i))^2 );
            
            if (d1-mrg<dis) && (dis<d1+mrg)
                HL(j,i) = t1;
                HL(i,j) = t1;
                
            elseif (d2-mrg<dis) && (dis<d2+mrg)
                HL(j,i) = t2;
                HL(i,j) = t2;
                
            elseif (d3-mrg<dis) && (dis<d3+mrg)
                HL(j,i) = t3;
                HL(i,j) = t3;
            end
        end
    end
end


% Generate HR
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

for j=1:N_Unit
    
    % Roll over atoms at the edges of AGNR to biuld Hamiltonian of ZCNT
    if ( (atom_y0(j) == 0) || (atom_y0(j) == Rows_Interval) )
        
        for i=1:N_Unit
            
            dis = sqrt( (atom_x0(j)-atom_xr(i))^2 + ((atom_y0(j)+Roll_Over)-atom_yr(i))^2 );
            
            if (d1-mrg<dis) && (dis<d1+mrg)
                HR(j,i) = t1;
                HR(i,j) = t1;
                
            elseif (d2-mrg<dis) && (dis<d2+mrg)
                HR(j,i) = t2;
                HR(i,j) = t2;
                
            elseif (d3-mrg<dis) && (dis<d3+mrg)
                HR(j,i) = t3;
                HR(i,j) = t3;
            end
        end
    end
end
fprintf('\t\t\t-> Done!\n\n');

disp=[];
p=1;
for ka=-pi:2*pi/500:pi
    HK=H0+HR*exp(1i*ka)+HL*exp(-1i*ka);
    ee=eig(HK);
    disp=[disp,ee];
    
    % Extract bandgap value
    for j=1:length(ee)
        if ee(j)>eps
            g(p)=ee(j);
            p=p+1;
        end
    end
end

%Band Gap
GAP=2*(min(g)-eps); if GAP< 1e-2 GAP=0; end

%Plot E-K diagram
figure();
plot(-pi:2*pi/500:pi,disp);
title('Band Structure');
xlabel('Ka');
ylabel('Energy(ev)');
string = sprintf('BG = %.2f eV, %d-ZCNT',GAP,n/2);
legend(string);
fprintf('\t\t-> done, Pristine Band Gap = %.2f eV\n\n',GAP);
fprintf('\n\n-> Done!\n\n');
fprintf('\n\n\n\t\t\tSimulation Is Complete');