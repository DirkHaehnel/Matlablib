function dy = CycloidElasticFun(t,y)

g = 9.81;
L = 1;
alpha = 0.1;

dy = zeros(2,1);
dy(1) = y(2);
dy(2) = -g*y(1)./(L + alpha*y(2));