function [qp, qs, pu, pd, su, sd, dpu, dpd, dsu, dsd] = HashWaveguideMode(n0,n,n1,d0,d,d1);

% function determines the in-plane wave vector components for the waveguide
% modes; it is assumed that n>max([n0 n1])

if n<max([n0 n1])
    disp('central layer is not waveguding');
    return
end
dq = 1e-4; % determines the finesse of searching; if waveguide modes are closer than this spacing one has a problem 
qq = max([n0(1) n1(end)])+dq/2:dq:max([n0 n n1])-dq/2;
w = sqrt(n^2-qq.^2);
rpu = Fresnel(w,[n n1],d1);
rpd = Fresnel(w,[n n0(end:-1:1)],d0(end:-1:1));
tst = unwrap(angle(rpu.*rpd.*exp(i*2*w*d)))/2/pi;
if tst(end)-tst(1)>0
    j1 = ceil(tst(1));
    j2 = floor(tst(end));
    dj = 1;
else
    j1 = floor(tst(1));
    j2 = ceil(tst(end));
    dj = -1;
end
qp = []; pu = []; pd = []; dpu = []; dpd = [];
for j = j1:dj:j2
    t = 1:length(qq);
    [t, t] = min(abs(tst-j));
    t = max([t-3, 1]):min([t+3, length(qq)]);
    qp = [qp; interp1(tst(t), qq(t), j, 'spline')];
    wp = sqrt(n^2-qp(end)^2);
    pu = [pu; interp1(qq(t), rpu(t), qp(end), 'spline')];
    pd = [pd; interp1(qq(t), rpd(t), qp(end), 'spline')];
    dpu = [dpu; interp1(w(t(1:end-1))+diff(w(t))/2, diff(rpu(t))./diff(w(t)), wp, 'spline')];    
    dpd = [dpd; interp1(w(t(1:end-1))+diff(w(t))/2, diff(rpd(t))./diff(w(t)), wp, 'spline')];
end

[rsu, rsu] = Fresnel(w,[n n1],d1);
[rsd, rsd] = Fresnel(w,[n n0(end:-1:1)],d0(end:-1:1));
tst = unwrap(angle(rsu.*rsd.*exp(i*2*w*d)))/2/pi;
if tst(end)-tst(1)>0
    j1 = ceil(tst(1));
    j2 = floor(tst(end));
    dj = 1;
else
    j1 = floor(tst(1));
    j2 = ceil(tst(end));
    dj = -1;
end
qs = []; su = []; sd = []; dsu = []; dsd = [];
for j = j1:dj:j2
    t = 1:length(qq);
    [t, t] = min(abs(tst-j));
    t = max([t-3, 1]):min([t+3, length(qq)]);
    qs = [qs; interp1(tst(t), qq(t), j, 'spline')];
    ws = sqrt(n^2-qs(end)^2);
    su = [su; interp1(qq(t), rsu(t), qs(end), 'spline')];
    sd = [sd; interp1(qq(t), rsd(t), qs(end), 'spline')];
    dsu = [dsu; interp1(w(t(1:end-1))+diff(w(t))/2, diff(rsu(t))./diff(w(t)), ws, 'spline')];    
    dsd = [dsd; interp1(w(t(1:end-1))+diff(w(t))/2, diff(rsd(t))./diff(w(t)), ws, 'spline')];
end

