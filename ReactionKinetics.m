close all

kp = 5;
km = 0.1;
K = kp/km;
a0 = 5;
b0 = 1; 
c0 = 1.2;

lam1 = -(b0+c0+K)/2;
lam2 = lam1 - sqrt(lam1^2 - b0*c0 + K*a0);
lam1 = lam1 + sqrt(lam1^2 - b0*c0 + K*a0);

t= 0:0.01:1;

x = lam1/lam2*exp(-(lam1-lam2)*km*t);
x = (lam1-lam2*x)./(1-x);

plot(t,a0-x,t,b0+x,t,c0+x)