 function [y, uscaled, start, decay, freq] = osc_plot_c(L, x, u)
% find the spatial eigenvalues
% roots of nu^4 - nu^2 + c == 0
% c is last element of u
c = u(end);
nu = roots([1 0 -1 0 c]);

% transfer this from [pi, 2 pi] to [0, L]
N  = length(x);
xL = linspace(0,L,N+1)';
xL = xL(1:end-1);

% oscillations frequency is imag(nu)
% oscillations decay with constant re(nu)
decay = abs(real(nu(1)));
freq  = abs(imag(nu(1)));

% where to start and end plot
l_bound = floor(length(xL)/50);
r_bound = length(xL) - floor(length(xL)/2);
y     = xL(l_bound:r_bound);

% scale solution by exp(decay) to recover oscillations
udata   = u(1:end-1);
uscaled = udata(l_bound:r_bound).*exp(decay*y);
umax    = max(uscaled(1:floor(end/2)));

% plot along with sine function of same scaling
figure;
plot(y,uscaled,y,umax.*sin(y*freq)); 

legend('rescaled solution','sine function')
title(strcat('scale to see oscillations, speed c =  ',num2str(c)))
start = l_bound
 
end
