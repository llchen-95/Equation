% The Crank Nicolson Method
% Example
% This CN solution uses theta=1/2
% File: CrankNicolsonP1.m
clc, clear, close
% Parameters
alfa=0.835;
dx=2;
dt=0.1;
nx=6; % or nx=[(10-0)/dx]+1 with L=10 cm
nt=11; % to compute 10 time steps k=[1,11]
% --- Constant Coefficients of the tridiagonal system
b = alfa/(2*dx^2); % Super diagonal: coefficients of u(i+1)
c = b; % Subdiagonal: coefficients of u(i-1)
a = 1/dt+b+c; % Main Diagonal: coefficients of u(i)
% Boundary conditions and Initial Conditions
Uo(1)=100; Uo(2:nx-1)=0; Uo(nx)=50;
Un(1)=100; Un(nx)=50;
% Store results for future use
UUU(1,:)=Uo;
% Loop over time
for k=2:nt

 for ii=1:nx-2
 if ii==1
 d(ii)=c*Uo(ii)+(1/dt-b-c)*Uo(ii+1)+b*Uo(ii+2)+c*Un(1);
 elseif ii==nx-2
 d(ii)=c*Uo(ii)+(1/dt-b-c)*Uo(ii+1)+b*Uo(ii+2)+b*Un(nx);
 else
 d(ii)=c*Uo(ii)+(1/dt-b-c)*Uo(ii+1)+b*Uo(ii+2);
 end
 end % note that d is row vector
 % Transform a, b, c constants in column vectors:
 bb=b*ones(nx-3,1);
 cc=bb;
 aa=a*ones(nx-2,1);
 % Use column vectors to construct tridiagonal matrices
 AA=diag(aa)+ diag(-bb,1)+ diag(-cc,-1);
 % Find the solution for interior nodes i=2,3,4,5
 UU=AA\d';
 % Build the whole solution as row vector
 Un=[Un(1),UU',Un(nx)];
 % Store results for future use
 UUU(k,:)=Un;
 % to start over
 Uo=Un;
end