function [xout, uout, c] = solveKdV_fourier_newton_interp(xin, uin, L, N)

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
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% setup

% compute finite difference matrices

[x, D]  = fourdif(N,1);
[x, D2] = fourdif(N,2);
[x, D3] = fourdif(N,3);
[x, D4] = fourdif(N,4);
[x, D5] = fourdif(N,5);
xout = x;

par.c = uin(end);
c = par.c;
u = interp1(xin,uin(1:end-1),xout);

% since we have no final grid point (periodic)
% will have some NaN values at the end of u
% replace these with 0
u(isnan(u)) = 0 ;

%% solve nonlinear problem using fsolve

% option to display output and use Jacobian
options=optimset('Display','iter','Jacobian','on','MaxIter',50);

% call solve
[uout,fval] = fsolve(@(u) integratedKdV_fourier(u,L,D,D2,D3,D4,D5,N,par),u,options);

%% plot results (optional)

% figure;
% plot(xout,uout); % plot solution
% title('Pulse on the full line');

% reappend c to the output vector
uout = [uout ; c];

