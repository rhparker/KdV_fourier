% generate initial data
% takes too long, so we will load it from a file
% N = 128;                 % number of finite diff gridpts
% iteration = 150;          % number of iterations of continuation code
% [x, uc] = solveKdV_fourier(N, iterations);

% instead we load initial data from file
% this contains 150 iterations with step 0.25 in c
%  x:  grid of x values, 0 to 50, 1000 grid points
%  uc: continuation data; each column of uc is of the form
%       [u; c], where u is the solution and c is the speed

%% load data and make configuration
load uc_50;

% method to use
config.method   = 'Fourier';
config.BC       = 'periodic';

% specify which index you want in uc
index       = 38;               % L = 50

% input data to use
u       = uc(:,index);

% ramp up to more grid points
% for increased accuracy
[xin, uin, c] = solveKdV_fourier_newton_interp(x, u, L, 1024);

% % if we don't do this, then use original values for input
% xin = x;
% uin = u;

% number of grid points is length of xin
N = length(xin);

% plot(xin,uin(1:end-1));

%% make plot of oscillations

% to find oscillations, we need half-wave
xhalf = xin(N/2+1:end);
uhalf = uin(N/2+1:end);

% plot oscillations
c = uin(end);
[y, uscaled, start, decay, freq] = osc_plot_c(L, xhalf, uhalf);

%% find the minima and maxima of this

% this does min/max using finite differences
% we can do this, since only finding min/max
udata   = uhalf(1:end-1);
% rescale xhalf to [0,L] like we did in the osc plot
xL = linspace(0,L,N/2+1)';
xL = xL(1:end-1);
uscaled = udata.*exp(decay*xL);
FDconfig.method   = 'fdiff';
FDconfig.BC       = 'Neumann';
N_diff = length(xL);
L_diff = xL(end) - xL(1);             
h      = L_diff/(N_diff - 1); 
D      = D_fdiff(N_diff, h, FDconfig.BC);
yhalf  = D*uscaled./uscaled - decay;

% threshold for a zero of this function
% this epsilon seems to work for this case
epsilon = 0.05;
% plot(xL,yhalf,'.');
% axis([0 L -epsilon, epsilon]);
zeros = find(abs(yhalf) < epsilon);

%% paste together half waves to make double pulses

% which zero to use
% 2, 3, 4, 5 are first four ones
i = 3;
z = zeros(i);

% construct the double wave
ud_left  = uin(z : N/2+z);
ud_right = flip(ud_left(1:end-1));
ud_right = ud_right(1:end-1);
ud_full  = [ ud_left; ud_right; c ];

% run double wave through Newton solver
ud_out  = solveKdV_fourier_newton(xin, ud_full, L);
figure;
plot(xin, ud_full(1:end-1), xin, ud_out(1:end-1));
legend('initial guess','Newton solver output')
title(strcat('double soliton, speed c =  ',num2str(c)))

% don't run this for now, since it takes forever
% and we have the output of eig saved

% [V,J] = eig_linear(xin, ud_out, config)