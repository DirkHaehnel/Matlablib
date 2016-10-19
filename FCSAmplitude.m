function amp = FCSAmplitude(t, x);
close all

p = simplex('rigler',[1e-4 1e-3 1e-5],[0 0 0],[],[],[],t,x);
[err, c, z] = rigler(p, t, x);
amp = c(1)/c(2);

