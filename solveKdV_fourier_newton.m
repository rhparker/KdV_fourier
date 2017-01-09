function [uout, c] = solveKdV_fourier_newton(xold, uold, L)

% - Solves 1D quadratic-cubic Swift-Hohenberg equation
% - Finds localised pulse in snaking region and computes its stability
% - Uses finite differences and sparse matrices
% - Requires optimization toolbox and external routine SH_rhs_finite_differences.m
%
% Copyright 2007, David JB Lloyd, Daniele Avitabile, Bjorn Sandstede, Alan R Champneys
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% D N
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% setup

% here use the Newton solver on the same domain as given
N = length(xold);

par.c = uold(end);
c = par.c;

u = uold(1:end-1);

%% compute Fourier differentiation matrices

[x, D]  = fourdif(N,1);
[x, D2] = fourdif(N,2);
[x, D3] = fourdif(N,3);
[x, D4] = fourdif(N,4);
[x, D5] = fourdif(N,5);

%% solve nonlinear problem using fsolve

% option to display output and use Jacobian
options=optimset('Display','iter','Jacobian','on','MaxIter',100);

% call solve
[uout,fval] = fsolve(@(u) integratedKdV_fourier(u,L,D,D2,D3,D4,D5,N,par),u,options);

%% plot results (optional)

% figure;
% plot(x,uout); % plot solution
% title('Pulse on the full line');

% reappend c to the output vector
uout = [uout ; c];

