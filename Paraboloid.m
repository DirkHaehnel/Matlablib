% Paraboloid

ng = 1.5;
nw = 1.333;
tc = [asin(nw/ng) [ceil(asin(nw/ng)):80]/180*pi];
tc = [0 asin(nw/ng)];

theta=(pi/10000:pi/5000:pi/2)';
z = (0:100)/100*5*pi;

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

        res(j,k) = real((vk-pk)/(vt-pt)+(vt*pk-pt*vk)*atan(sqrt((vt-pt)/pt))/sqrt(pt)/(vt-pt)^(3/2));        
    end
end

plot(tc/pi*180,res(1,:),[tc(1) tc(1)]/pi*180,[0 0.35],':')
xlabel('angle [°]');
ylabel('collection efficiency');
figure
ind = [1 2 3 4 5 7 10 19];
h=plot(z/2/pi,res(:,ind)./(ones(size(res,1),1)*max(res(:,ind))))
colorize(h)
legend(num2str(tc(ind)'/pi*180))
xlabel('distance [\lambda]');
ylabel('rel. collection efficiency');

