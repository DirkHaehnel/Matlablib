% azimut angle phi 
phi = 0:1/180*pi:360/180*pi; phi = phi';
% polar angle  thet
thet = -90/180*pi:0.01/180*pi:90/180*pi;

% refractive indexes
n1 = 1.33; n2 =1.5;

% surface distance in vacuum wavelength

for j = 1:61,

z0 = (j-1)/40*pi*2;

p = 1;
e1 = n1^2; e2 = n2^2;

% orthogonal dipol
% upper half space
q = n1*sin(thet); w1 = sqrt(n1^2-q.^2); w2 = sqrt(n2^2-q.^2);
rp = (w1*e2-w2*e1)./(w1*e2+w2*e1);
s = n1^3/(e1^2)*q.^2.*(abs(1+rp.*exp(2*i*w1*z0))).^2;

% orthogonal dipol
% lower half space
q = n2*sin(thet); w1 = sqrt(n1^2-q.^2); w2 = sqrt(n2^2-q.^2);
tp = 2*n1*n2*w1./(w1*e2+w2*e1);
t = n2*n1^2*w2.^2./(e1^2*(abs(w1)).^2).*q.^2.*(abs(tp)).^2.*exp(-2*imag(w1)*z0);

% parallel dipol 
% upper half space
q = n1*sin(thet); w1 = sqrt(n1^2-q.^2); w2 = sqrt(n2^2-q.^2);
rp = (w1*e2-w2*e1)./(w1*e2+w2*e1);
rs = (w1-w2)./(w1+w2);
u = ones(length(phi),length(thet))*n1.*(((cos(phi)).^2*(((abs(1-rp.*exp(2*i*w1*z0))).^2).*w1.^2/n1^2))+(sin(phi)).^2*(abs(1+rs.*exp(2*i*w1*z0))).^2);

% parallel dipol 
% lower half space
q = n2*sin(thet); w1 = sqrt(n1^2-q.^2); w2 = sqrt(n2^2-q.^2);
tp = 2*n1*n2*w1./(w1*e2+w2*e1);
ts = 2*w1./(w1+w2);
v = ones(length(phi),1)*(n2*w2.^2./(abs(w1)).^2.*exp(-2*imag(w1)*z0)).*(((cos(phi)).^2*((abs(w1)).^2/n1^2.*(abs(tp)).^2)+(sin(phi)).^2*(abs(ts)).^2));

% averaging the azimut angles:
um = mean(u); vm = mean(v);

% isotropic orientation:   
upp =2*um;
low =2*vm;


% factor sin(theta) 
upper=upp.*abs(cos(thet-pi/2));
lower=low.*abs(cos(thet-pi/2));

hh(j)=sum(upper)/(sum(upper)+sum(lower));
end;

k1=polar(thet-pi/2,low,'k');
set(k1,'LineWidth',2,'Color','red');
hold on;
k2=polar(thet+pi/2,upp,'k');
set(k2,'Linewidth',2,'Color','red');
hold off;



   
        
       

