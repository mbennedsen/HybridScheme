function X = simBSS(a,Lfct,n,T,gamm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate BSS process observed on equidistant grid using the Hybrid Scheme of Bennedsen, Lunde, and Pakkanen (2017). 
%
% INPUT
% a          : Roughness parameter; alpha \in (-0.5,0.5).
% Lfct       : Function handle of the slowly varying function to include (see below for details).
% n          : Number og observations to simulate.
% T          : Terminal time. I.e. the process is simulated on the grid dt, 2*dt, ... , n*dt = T, where dt = T/n.
% gamm       : Parameter used in hybrid scheme to set "burn in" period. 
%
% OUTPUT
% X          : (n x 1) vector of simulated observations of the BSS process.
%
% NOTE
% In this version, the value kappa=3 is hard-coded. That is, the first
% three terms in the Riemann-sum approximation to the BSS process have been
% approximated using the Hybrid method.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLOWLY VARYING FUNCTION
% The kernel function of the simulated BSS process is g(x) = x^a*Lfct(x),
% where Lfct(x) is an input function handle containing a slowly varying function. 
%
%
% EXAMPLES OF SLOWLY VARYING FUNCTIONS TO USE
%
% Gamma-BSS Process (Example 2.3):
% Lfct = @(x)( exp( -lambda*x ) ), with lambda>0.
%
% Power-BSS Process (Example 2.4):
% Lfct = @(x)( (1+x).^(beta-a) ), with beta<-0.5. (Note that this process
% has the long memory property for beta \in (-1,-0.5), see Example 2.4 in
% the paper.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (c) Mikkel Bennedsen (2021)
%
% This code can be used, distributed, and changed freely. Please cite Bennedsen,
% Lunde, and Pakkanen (2017): "Hybrid scheme for Brownian semistationary processes", Finance and Stochastics (2017), 21, 931-965.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks
if nargin < 4
    T = 1;
end

if nargin < 5
    gamm = 0.5;
end


%% Init
dt = T/n;
N = floor(n^(1+gamm));

%% Stuff for simulating BSS
b = ( ((1:N).^(a+1) - ((1:N)-1).^(a+1))/(a+1)).^(1/a) - (1:N);
t1 = (1:N)*dt;
t2 = t1+b*dt;
g   = (t2.^a).*Lfct(t2);%exp(-l*t2);

%% Construct appropriate (1+kappa)x(1+kappa) matrix. Here, kappa=3.
covM = nan(4,4);

% Diagonal
covM(1,1) = dt;
covM(2,2) = dt^(2*a+1)/(2*a+1);
covM(3,3) = dt^(2*a+1)*(2^(2*a+1) - 1)/(2*a+1);
covM(4,4) = dt^(2*a+1)*(3^(2*a+1) - 2^(2*a+1))/(2*a+1);

% Row 1
covM(1,2) = dt^(a+1)/(a+1);                     covM(2,1) = covM(1,2);
covM(1,3) = dt^(a+1)*(2^(a+1) - 1)/(a+1);       covM(3,1) = covM(1,3);
covM(1,4) = dt^(a+1)*(3^(a+1) - 2^(a+1))/(a+1); covM(4,1) = covM(1,4);

% Row 2
covM(3,2) = dt^(2*a+1)*( (2*1)^(a+1)*hypergeom([1,2*(a+1)],a+2,-1))/(a+1);          covM(2,3) = covM(3,2);
covM(4,2) = dt^(2*a+1)*( (3*1)^(a+1)*hypergeom([1,2*(a+1)],a+2,-1/2))/(3-1)/(a+1);  covM(2,4) = covM(4,2);

% Row 3
covM(4,3) = dt^(2*a+1)*( (3*2)^(a+1)*hypergeom([1,2*(a+1)],a+2,-2) - (2*1)^(a+1)*hypergeom([1,2*(a+1)],a+2,-1))/(3-2)/(a+1);   covM(3,4) = covM(4,3);

FF2 = chol(covM);

%% Do simulation
We = FF2'*randn(4,n+N);
dW1 = We(1,:);

%% Hybrid convolution step
X = conv([0,0,0,g(4:end)],dW1);
X = X(N+1:N+n) + We(2,N+1:N+n)*Lfct(t2(1)) + We(3,N:N+n-1)*Lfct(t2(2)) + We(4,N-1:N+n-2)*Lfct(t2(3));
X = X';

