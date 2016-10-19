function rad = Diff2Rad(dc,tt,vv)

if length(tt)<length(dc)
    tt = tt*ones(size(dc));
end
tt(tt<100) = tt(tt<100) + 273.15;
T = 273.15 + (0:10:100);   
BoltzmannConstant = 1.380650324e-23; % NIST
if nargin<3 || isempty(vv)
    visc = [1.7920 1.3080 1.0050 0.8007 0.6560 0.5494 0.4688 0.4061 0.3565 0.3165 0.2838];
    vv = interp1(T,visc,tt,'cubic');
end
rad = tt*BoltzmannConstant/6/pi./vv./dc/1e-7*1e10; % hydrodyn. radius in Angstrom
