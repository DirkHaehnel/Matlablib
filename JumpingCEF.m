ng = 1.5;
nw = 1.;
tc = [0 asin(nw/ng)];

theta=(pi/10000:pi/5000:pi/2)';
z = (0:100)/100*pi;

st = sin(theta);
ct = cos(theta);

for j=1:length(z)
    for k=1:length(tc)
        [v1,pc1,ps1] = DipoleL(theta,z(j),ng,nw,nw,[],10,[]);
        [v2,pc2,ps2] = DipoleL(theta,10-z(j),nw,nw,ng,[],10,[]);
        
        sk = sin(theta(theta>tc(k)));
        ck = cos(theta(theta>tc(k)));
        
        vk = sum(sk.*abs(v1(theta>tc(k))).^2);
        pk = sum(sk.*(abs(pc1(theta>tc(k))).^2+abs(ps1(theta>tc(k))).^2))/2;
        
        vt = sum(st.*(abs(v1).^2+abs(v2).^2));
        pt = sum(st.*(abs(pc1).^2+abs(ps1).^2+abs(pc2).^2+abs(ps2).^2))/2;

        rva(j,k) = vk/vt;
        rpa(j,k) = pk/pt;
        resa(j,k) = real((vk-pk)/(vt-pt)+(vt*pk-pt*vk)*atan(sqrt((vt-pt)/pt))/sqrt(pt)/(vt-pt)^(3/2));        
        
        [v1,pc1,ps1] = DipoleL(theta,5*pi-z(j),ng,ng,nw,[],5*pi,[]);
        [v2,pc2,ps2] = DipoleL(theta,z(j),nw,ng,ng,[],10,[]);
        
        sk = sin(theta(theta>tc(k)));
        ck = cos(theta(theta>tc(k)));
        
        vk = sum(sk.*abs(v1(theta>tc(k))).^2);
        pk = sum(sk.*(abs(pc1(theta>tc(k))).^2+abs(ps1(theta>tc(k))).^2))/2;
        
        vt = sum(st.*(abs(v1).^2+abs(v2).^2));
        pt = sum(st.*(abs(pc1).^2+abs(ps1).^2+abs(pc2).^2+abs(ps2).^2))/2;

        rvb(j,k) = vk/vt;
        rpb(j,k) = pk/pt;
        resb(j,k) = real((vk-pk)/(vt-pt)+(vt*pk-pt*vk)*atan(sqrt((vt-pt)/pt))/sqrt(pt)/(vt-pt)^(3/2));        
    end
end

h=plot([-fliplr(z),z]/2/pi,[flipud(resb); resa]./(ones(2*size(resa,1),1)*max([resa;resb])))
colorize(h)
ax = axis;
hold on
plot([0 0],[0 ax(4)],':g')
hold off
legend(num2str(tc'/pi*180))
xlabel('distance [\lambda]');
ylabel('rel. collection efficiency');

figure
h=plot([-fliplr(z),z]/2/pi,[flipud(rvb); rva]./(ones(2*size(rva,1),1)*max([resa;resb])))
ax = axis;
hold on
plot([0 0],[0 ax(4)],':g')
hold off
colorize(h)
legend(num2str(tc'/pi*180))
xlabel('distance [\lambda]');
ylabel('rel. collection efficiency for vertical dipole only');

figure
h=plot([-fliplr(z),z]/2/pi,[flipud(rpb); rpa]./(ones(2*size(rpa,1),1)*max([resa;resb])))
ax = axis;
hold on
plot([0 0],[0 ax(4)],':g')
hold off
colorize(h)
legend(num2str(tc'/pi*180))
xlabel('distance [\lambda]');
ylabel('rel. collection efficiency for parallel dipole only');
