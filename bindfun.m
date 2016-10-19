function [err, y, tmp] = bindfun(p, conc, int);

conc = conc(:)';
int = int(:)';
n = (0:100)';

tmp = (conc + p(1) + p(2))/2; 
tmp = tmp - sqrt(tmp.^2 - p(1)*conc);

y = tmp.*exp(-(tmp/p(3)).^p(4));
%y = tmp./(1+(tmp/p(3)).^p(4));

%y = abs(int/[ones(size(conc)); conc-tmp; y])*[ones(size(conc)); conc-tmp; y];
y = (int/[ones(size(conc)); y])*[ones(size(conc)); y];

plot(conc, int, 'o', conc, y); drawnow;
%err = sum((int - y).^2.*abs(y));
err = abs(sum((int - y).^2));
disp([err p'])