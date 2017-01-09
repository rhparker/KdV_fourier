% eigenvalues of linearization about stationary wave
% uses eigs instead of eig for speed
function [lambda, V, J] = eigs_linear(x, u, L, num, center)
    par.c = u(end);
    udata = u(1:end-1);
    N = length(udata);
    [~, D]  = fourdif(N,1);
    [~, D2] = fourdif(N,2);
    [~, D3] = fourdif(N,3);
    [~, D4] = fourdif(N,4);
    [~, D5] = fourdif(N,5);
    [F, J] = KdV_fourier(udata,L,D,D2,D3,D4,D5,N,par);
    opts.tol = 10^(-10); 
    [V, lambda_D] = eigs(J, num, center, opts);
    lambda = diag(lambda_D);
end
