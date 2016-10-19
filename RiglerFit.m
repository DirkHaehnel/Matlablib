function [dc a F triplet c err z] = RiglerFit(data, p, m, bounds)


if (nargin < 3)  || (isempty(m))
    m = 0;
end

if (nargin < 2)
    p = [4e-3 3.3 1/1e-5];
else
    if (m == 1)
        p(3) = p(3) * 2e6 * pi();
    end
    p(3+m:end) = 1./p(3+m:end);
end

if (nargin < 4) || isempty(bounds)
	xmin = zeros(1, length(p));
    xmax = [];
else
    xmin =  bounds(1,:);
    xmax =  bounds(2,:);
end

[x, dx, steps] = Simplex('Rigler_Cos',p,xmin,xmax,[],[],data.t,data.y,m,1);
[err, c, z] = Rigler_Cos(x,data.t,data.y,m,1);

ylabel('correlation / (\ith\rm\nu)^2/s^2'); 
xlabel('time (s)')

s = [];
s{1} = ['\tau_D = ' num2str(x(1)*1e3,3) ' ms'];
s{2} = ['a = ' num2str(x(2),3)];

dc = x(1)*1e3;
a = x(2);

if (m > 0)
    F = x(3:2+m) / 2e6 / pi();
    s{end+1} = ['F = ' num2str(x(3),4) 'kHz'];
else
    F = [];
end

triplet = 1./x(3+m:end)*1e6;

if (length(p)-2-m > 0)
    tmp = ['\tau_t = (' num2str(1/x(3+m)*1e6,4)];
    for k=4+m:length(p)
        tmp = [tmp ',' num2str(1/x(k)*1e6,4)];
    end
    tmp = [tmp  ') \mus'];
    s{end+1} = tmp;
end

text(0.6,0.95,s,'VerticalAlignment','Top','units','normal')