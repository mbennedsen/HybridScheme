%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Illustration of how to use the simBSS function, which simulates a BSS process using the Hybrid Scheme of Bennedsen, Lunde, and Pakkanen (2017).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (c) Mikkel Bennedsen (2021)
%
% This code can be used, distributed, and changed freely. Please cite Bennedsen,
% Lunde, and Pakkanen (2017): "Hybrid scheme for Brownian semistationary processes", Finance and Stochastics (2017), 21, 931-965.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;
%% Init
a      = -0.25;    % Roughness index
lambda = 0.1;      % Memory parameter for Gamma-BSS process (lambda>0).
beta   = -0.75;    % Memory parameter for Power-BSS process (beta<-0.5).

T = 100;           % Terminal time
n = 1000;          % Number of observations to simulate

gamm = 0.50;

%% L fct
Lfct1 = @(x)( exp(-lambda*x) );  % Slowly varying function for Gamma-BSS
Lfct2 = @(x)( (1+x).^(beta-a) ); % Slowly varying function for Power-BSS

%% Stuff
dt = T/n;
t = (dt:dt:T)';

%% Sim processes
X1 = simBSS(a,Lfct1,n,T,gamm); % Gamma-BSS
X2 = simBSS(a,Lfct2,n,T,gamm); % Power-BSS

%% Plot
figure;
plot(t,X1,'LineWidth',1.5), hold on
plot(t,X2,'LineWidth',1.5), hold on

legend('Gamma-BSS','Power-BSS','Location','Best');