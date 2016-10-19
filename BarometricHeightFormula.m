NatConst
mo2 = 2*16e-3;
mn2 = 2*14e-3;
xo2 = 0.22;
xn2 = 0.78;
h = 2e3;
t = IcePoint + 25;
p0 = 100e3;

p = p0*(xo2*exp(-mo2*AccelerationDueToGravity*h/MolarGasConstant/t) + ...
    xn2*exp(-mn2*AccelerationDueToGravity*h/MolarGasConstant/t))
