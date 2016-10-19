function y = doubfit(kap, data, bb);

mx = size(bb);
mx = mx(1);
tmp = exp((0:(mx-1))*log(1-kap))*kap^2;
tmp = tmp*bb;
c = data/tmp;
y =sum((data - c*tmp).^2);

