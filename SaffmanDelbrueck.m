function [err, c, z] = SaffmanDelbrueck(p,R,dif,etaf,T,flag)

R = R(:);
kB =  1.38065e-023;
gamma = 0.5772156649;
z = kB*T/(4*pi*p)*(log(p/etaf./R)-gamma);

if length(dif)==1
    err = dif*z;
else
    dif = dif(:);
    if nargin<6 || isempty(flag) % keep amplitude fixed or not
    	c = z\dif;
        z = c*z;
    end
    
    plot(R,dif,'o',R,z); drawnow
    
    err = sum((dif-z).^2./z);
end

return

R = 1e-9*[0.4 0.5 1 1.56 1.8 2.1 2.4 2.8 3.6]';
Ri = linspace(min(R),max(R),100);

dif = 1e-12*[11.575 10.30375 9.272 8.985 8.49 8.46556 7.955 7.5125 7.22667]';
der = 1e-12*[0.56789 0.60738 0.53364 0.8338 0.46476 0.7488 0.5252 0.83994 0.39272];


load Viscosity
tt = 21;
etaf = interp1(T,visc(1,:),273.15+tt)*1e-3;

p = Simplex('SaffmanDelbrueck',etaf*5e-7,[],[],[],[],R,dif,etaf,273.15+tt);
[err, c] = SaffmanDelbrueck(p,R,dif,etaf,273.15+tt);
y1i = SaffmanDelbrueck(p,Ri,c,etaf,273.15+tt);

c = (1./R)\dif;
y2i = c./Ri;

plot(R*1e9,dif*1e12,'s',Ri*1e9,y1i*1e12,Ri*1e9,y2i*1e12,'MarkerSize',5,'MarkerFaceColor','r')
for j=1:length(R)
    line(1e9*R(j)*[1 1],1e12*(dif(j)+[-der(j) der(j)]),'Color','r')
end
xlabel('radius (nm)'); 
ylabel('diffusion coefficient (\mum^2/s)')

