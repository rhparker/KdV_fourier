% piece together two half-waves to get a full wave

function [xfull, ufull] = full_wave(xin, uin)
% remove speed c from the end
udata = uin(1:end-1);
c = uin(end);
% piece together two halves to get whole wave
uleft = flipud(udata);
ufull = [ uleft(1:end-1) ; udata; c ];

% piece together domain
xleft = -flipud(xin);
xfull = [ xleft(1:end-1); xin ];

end