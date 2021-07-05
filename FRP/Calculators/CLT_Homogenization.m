% =============================================== %
% Composite homogenization using Composite Lamination Theory (CLT)
% Author: Paul Davidson, paul.davidson@uta.edu
% =============================================== %
clear all;

% LAMINA INFORMATION

N=1;            % Number of repeating lamina
M=1;            % Number of material systems

% Layers 25/50/25, 33/67/0
sym=0;          % Is there symmetry?
t=0.25*[1,2,1];    % Stacking order - thickness (in mm)
a=[0,90,0];    % Stacking order - Angle (in degrees)
%a=[60,-60,60,0,0,-60,60,0,-60,60,-60,0];     % Stacking order - Angle (in degrees)
e=[1,2,1];    % Stacking order - Material

 
a=a*pi/180; % angle deg to radian
T=sum(t)*N*2^(sym);       % Total thickness

% MATERIAL INFORMATION
E1=[169,45];
E2=[9,11];
v12=[0.31,0.29];
v21=[v12(1)*E2(1)/E1(1),v12(2)*E2(2)/E1(2)];
G12=[6.5,4.5];

% Q 
i=1;
for j=1:1:length(t)*N
    
    km=e(i);
    L(j).G=[1 0 0;0 1 0;0 0 1];
    L(j).S=[ 1/E1(km) -v21(km)/E2(km) 0;
            -v12(km)/E1(km) 1/E2(km) 0;
            0 0 1/G12(km)];
    L(j).Q=[ (E1(km)/(1-v12(km)*v21(km))) (v12(km)*E2(km)/(1-v12(km)*v21(km))) 0;
        (v12(km)*E2(km)/(1-v12(km)*v21(km))) (E2(km)/(1-v12(km)*v21(km))) 0;
        0 0 G12(km)];
    L(j).T1=[ cos(a(i))^2 sin(a(i))^2 2*cos(a(i))*sin(a(i));
        sin(a(i))^2 cos(a(i))^2 -2*cos(a(i))*sin(a(i));
        -cos(a(i))*sin(a(i)) cos(a(i))*sin(a(i)) (cos(a(i))^2 - sin(a(i))^2)];
    L(j).T2=[ cos(a(i))^2 sin(a(i))^2 cos(a(i))*sin(a(i));
        sin(a(i))^2 cos(a(i))^2 -cos(a(i))*sin(a(i));
       -2*cos(a(i))*sin(a(i)) 2*cos(a(i))*sin(a(i)) (cos(a(i))^2 - sin(a(i))^2)];    
    L(j).Qb=(L(i).T1) \ L(i).Q*(L(i).T2);
    L(j).t=t(i);
    L(j).Sb=inv(L(j).Qb);
    % Repeat count
    if i>=length(t)
        i=1;
    else
        i=i+1;
    end
    
end

% IMPLEMENTATION OF SYM
if sym==1
    for j=1:1:length(t)*N
        L(length(t)*N+j)=L(length(t)*N-j+1);
    end
end

% CALCULATION OF H
    he=0;
    for k=1:1:length(t)*N*2^(sym)
    L(k).H=[he-T/2 he+L(k).t-T/2];
    he=he+L(k).t;
    end
    
    
% LAMINATE MATRIX
A=0;
B=0;
D=0;
for i=1:1:length(L)
    A=A+L(i).Qb*(L(i).H(2)-L(i).H(1));
    B=B+(0.5*L(i).Qb*(L(i).H(2)^2-L(i).H(1)^2));
    D=D+((1/3)*L(i).Qb*(L(i).H(2)^3-L(i).H(1)^3));    
end

% ZERO OUT SMALL NUMBERS
for i=1:1:3
    for j=1:1:3
        if A(i,j)<1e-6 && A(i,j)>-1e-6
           A(i,j)=0;
        end
        if B(i,j)<1e-6 && B(i,j)>-1e-6
           B(i,j)=0;
        end
        if D(i,j)<1e-6 && D(i,j)>-1e-6
           D(i,j)=0;
        end
    end
end

disp('`````````````````````````````````````````````')
disp('ABD Matrices')
A
B
D
disp('`````````````````````````````````````````````')
B_str=-inv(A)*B;
C_str=B*inv(A);
D_str=D-B*inv(A)*B;

disp('a b c d Matrices')
a=inv(A)-(B_str*inv(D_str))*C_str
b=B_str*inv(D_str)
c=-inv(D_str)*C_str
d=inv(D_str)
% Laminate constants
disp('--------------------')
disp('Homogenized Laminate Properties')
As=T*inv(A);
Ex=1/As(1,1)
Ey=1/As(2,2)
muxy=-As(1,2)/As(1,1)
Gxy=1/As(3,3)
disp('--------------------')
Ex=(A(1,1)-A(1,2)^2/A(2,2))*1/T
Ey=(A(2,2)-A(1,2)^2/A(1,1))*1/T
muxy=A(1,2)/A(2,2)
Gxy=A(3,3)/T