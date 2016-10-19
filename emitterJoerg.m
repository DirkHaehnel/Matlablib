% polar angle  theta
theta = pi/1000:pi/500:pi/2;

% refractive indexes
n1 = 1.33; n2 =1.5;

z0 = (0:0.2:60)/40*pi*2;

dpw = zeros(length(z0),1); upw = dpw; dvw = dpw; uvw = dpw;
dpg = dpw; upg = dpg; dvg = dpg; uvg = dpw;
for j=1:length(z0)
    [v,pc,ps] = DipoleL(theta,z0(j),n2,n1,n1,[],0,[]);
    dpw(j) = sin(theta)*(abs(pc).^2+abs(ps).^2)/2;
    dvw(j) = sin(theta)*(abs(v).^2);
    
    [v,pc,ps] = DipoleL(theta,max(z0)-z0(j),n1,n1,n2,[],max(z0),[]);
    upw(j) = sin(theta)*(abs(pc).^2+abs(ps).^2)/2;
    uvw(j) = sin(theta)*(abs(v).^2);    

    [v,pc,ps] = DipoleL(theta,max(z0)-z0(j),n2,n2,n1,[],max(z0),[]);
    dpg(j) = sin(theta)*(abs(pc).^2+abs(ps).^2)/2;
    dvg(j) = sin(theta)*(abs(v).^2);
    
    [v,pc,ps] = DipoleL(theta,z0(j),n1,n2,n2,[],max(z0),[]);
    upg(j) = sin(theta)*(abs(pc).^2+abs(ps).^2)/2;
    uvg(j) = sin(theta)*(abs(v).^2);    

end

plot([-fliplr(z0) z0],[flipud(dpg./(dpg+upg)); dpw./(dpw+upw)],[-fliplr(z0) z0],[flipud(dvg./(dvg+uvg)); dvw./(dvw+uvw)])

plot([-fliplr(z0) z0],[flipud(upg./(dpg+upg)); upw./(dpw+upw)],[-fliplr(z0) z0],[flipud(uvg./(dvg+uvg)); uvw./(dvw+uvw)],zeros(101,1),0.1+(0:100)/100*0.55,':')
axis tight
xlabel('distance (nm)');
ylabel('rel. emission into water');