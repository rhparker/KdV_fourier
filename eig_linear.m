% eigenvalues of linearization about stationary wave
function [lambda, V, J] = eig_linear(x, u, L)
    par.c = u(end);
    udata = u(1:end-1);
    N = length(udata); 
    [~, D]  = fourdif(N,1);
    [~, D2] = fourdif(N,2);
    [~, D3] = fourdif(N,3);
    [~, D4] = fourdif(N,4);
    [~, D5] = fourdif(N,5);
    [F, J] = KdV_fourier(udata,L,D,D2,D3,D4,D5,N,par);
    [V, LD] = eig(full(J));
    lambda = diag(LD);
end
