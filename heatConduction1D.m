%% Homework #2
%% Solution to 1-D heat transfer using 
%% Tri-diagonal algorithm
%% Anurag Sandeep K. 
%% EQ. -d/dt(kAdT/dx)+beta*P*(T-Tinf)=0
%% Dependencies - Tridiagonal.m
% k - thermal conductivity [w/mC]
% area - Area of cross-section
% beta
% P - Perimeter
% A - co-efficient matrix
% B - solution vector
% T1 - T0-Tinf = 300
% N = number of grid points=(ITMAX)-2
% deltax - grid length
clc
clear all
%-------------------------------------
% INPUT PARAMETERS
%-------------------------------------
k=50; 
D=0.02;
L=1.0;
area=pi*D^2/4;
beta=100;
P=pi*D;
T1=300;
ITMAX=[6 11 21 41 81];
A=[];B=[];T=[];
%-------------------------------------
% MAIN LOOP TO COMPUTE A & B MATRICES
% and call TDMA Algorithm
%-------------------------------------
for i=1:length(ITMAX)
    N(i)=ITMAX(i)-2;
    deltax(i)=(L-0)/(N(i)+1);

    % Equation #1 - 2nd node from left
    A(1,2)=k*area/(deltax(i)/2) + k*area/deltax(i) + beta*P*deltax(i);
    A(1,3)=-k*area/deltax(i);
    B(1)=k*area*T1/(deltax(i)/2);

    % Equation #2 - Internal nodes
    for IT=3:N(i)
        A(IT-1,1)=-k*area/deltax(i);
        A(IT-1,2)=2*k*area/(deltax(i)) + beta*P*deltax(i);
        A(IT-1,3)=-k*area/deltax(i);
        B(IT-1)=0;   
    end
    
    % Equation #N-1 - (n-1)th node
    A(N(i),1)=-k*area/deltax(i);
    A(N(i),2)=2*k*area/(deltax(i)) + k*area/deltax(i) + beta*P*deltax(i);
    A(N(i),3)=-k*area/deltax(i)/2;
    B(N(i))=0;

    % Equation #N - for the right most boundary node
    A(N(i)+1,1)=-k*area/(deltax(i)/2) ;
    A(N(i)+1,2)=k*area/(deltax(i)/2)+ beta*area;
    B(N(i)+1)=0;

% Calling Tri-diagonal algorithm
b=Tridiagonal(N(i),A,B);
b=b+20; 
clear A B

% Plot temperature
plot(b,'LineWidth',2); hold on
xlabel('$ITMAX$','Interpreter','latex');
ylabel('$T_{(2:ITMAX)}$','Interpreter','latex');
legend('ITMAX=6','ITMAX=11','ITMAX=21','ITMAX=41','ITMAX=81');
title('Temp. distribution using Tri-diagonal Matrix Algorithm')
end






