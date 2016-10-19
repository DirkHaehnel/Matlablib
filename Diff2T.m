function tt = Diff2T(dc,rad)

T = 273.15 + (0:10:100);   
visc = [1.7920 1.3080 1.0050 0.8007 0.6560 0.5494 0.4688 0.4061 0.3565 0.3165 0.2838];
BoltzmannConstant = 1.380650324e-23; % NIST
dd = T*BoltzmannConstant/6/pi./visc./rad/1e-7*1e10; 
tt = interp1(dd,T-273.15,dc,'cubic');
