function [err, c, y] = KinFun(p, conc, t, int, bld);

% Resact to cGMP kinetics
% parameter values: 
% p(1) - max receptors
% p(2) - binding constant
% p(3) - guanylate cyclase activity
% p(4) - guanylate cyclase inhibition
% p(5) - inhibition coefficient

conc = conc(:)';
t = t(:);
y = zeros(length(t),length(conc));
for j = 1:length(conc)
   [tmp, tmp] = ode15s(@kinetics, t, [0 0 0], [], p, conc(j));
   if size(tmp,1)==size(y,1)
       y(:,j) = tmp(:,3); 
   else
       y(:,j) = 0; 
   end
end

col = ones(size(y(:)));
err = 0;
if nargin>3 & ~isempty(int) & sum(y(:))>0
    if ~size(t,1)==size(int,1) int=int'; end
    ind = ~isnan(int);
    c = y(ind)\int(ind);
    % c = 1;
    z = c*y;
    err = sum(sum(abs(z(ind)-int(ind)).^2./(abs(z(ind))+(abs(z(ind))==0))));
    if nargin>4 & ~isempty(bld)
        plot(t,int,'o'); hold on; plot(t,z); hold off; drawnow    
    end
else
    err = inf;
end

% --------------------------------------------------------------------------

function dy = kinetics(t,y,p,conc);

dy = zeros(3,1);
dy(1) = p(2)*(conc-y(1))*(p(1)-y(1));
if length(p)==4
    dy(2) = p(2)*(conc-y(1))*(p(1)-y(1)) - p(4)*y(2);
else
    dy(2) = p(2)*(conc-y(1))*(p(1)-y(1)) - p(4)*y(2)*y(3)^p(5);
end
dy(3) = p(3)*y(2);



