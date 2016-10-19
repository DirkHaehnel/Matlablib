function [fr, ff, fz] = CylField(nv,rv,m,w,coef,rr,phi)

fr = []; ff = []; fz = [];
ind = rr>=rv(end);
if sum(ind)>0
    [br, bf, bz] = CylN(m,sqrt(nv(end)^2-w.^2),w,rr(ind),'h');
    [ar, af] = CylM(m,sqrt(nv(end)^2-w.^2),rr(ind),'h');
    fr = coef(1)*ar+coef(2)*br;
    ff = coef(1)*af+coef(2)*bf;
    fz = coef(2)*bz;
end
for j=length(rv):-1:2
    ind = rr>=rv(j-1) & rr<rv(j);
    if sum(ind)>0
        [br, bf, bz] = CylN(m,sqrt(nv(j)^2-w.^2),w,rr(ind));
        [ar, af] = CylM(m,sqrt(nv(j)^2-w.^2),rr(ind));
        [dr, df, dz] = CylN(m,sqrt(nv(j)^2-w.^2),w,rr(ind),'h');
        [cr, cf] = CylM(m,sqrt(nv(j)^2-w.^2),rr(ind),'h');
        fr = [coef(4*(length(rv)-j)+3)*ar+coef(4*(length(rv)-j+1))*br+coef(4*(length(rv)-j+1)+1)*cr+coef(4*(length(rv)-j+1)+2)*dr fr];
        ff = [coef(4*(length(rv)-j)+3)*af+coef(4*(length(rv)-j+1))*bf+coef(4*(length(rv)-j+1)+1)*cf+coef(4*(length(rv)-j+1)+2)*df ff];
        fz = [coef(4*(length(rv)-j+1))*bz+coef(4*(length(rv)-j+1)+2)*dz fz];
    end
end
ind = rr<=rv(1);
if sum(ind)>0
    [br, bf, bz] = CylN(m,sqrt(nv(1)^2-w.^2),w,rr(ind));
    [ar, af] = CylM(m,sqrt(nv(1)^2-w.^2),rr(ind));
    fr = [coef(end-1)*ar+coef(end)*br fr];
    ff = [coef(end-1)*af+coef(end)*bf ff];
    fz = [coef(end)*bz fz];
end

if nargin>6 && ~isempty(phi)
    phi = phi(:);
    fr = exp(i*m*phi)*fr;
    ff = exp(i*m*phi)*ff;
    fz = exp(i*m*phi)*fz;
end

% subplot(131); surf(cos(phi)*rr,sin(phi)*rr,abs(fr).^2); shading interp; view([-36 10]); grid off; light
% subplot(132); surf(cos(phi)*rr,sin(phi)*rr,abs(ff).^2); shading interp; view([-36 10]); grid off; light
% subplot(133); surf(cos(phi)*rr,sin(phi)*rr,abs(fz).^2); shading interp; view([-36 10]); grid off; light
% colormap hot
% 
% subplot(131); pcolor(cos(phi)*rr,sin(phi)*rr,abs(fr).^2); axis image; shading interp; colorbar('h'); 
% subplot(132); pcolor(cos(phi)*rr,sin(phi)*rr,abs(ff).^2); axis image; shading interp; colorbar('h'); 
% subplot(133); pcolor(cos(phi)*rr,sin(phi)*rr,abs(fz).^2); axis image; shading interp; colorbar('h'); 
% colormap hot
% 
% plot(rr/2/pi,imag(fr(1,:))+real(fr(1,:)),rr/2/pi,imag(ff(1,:))+real(ff(1,:)),rr/2/pi,imag(fz(1,:))+real(fz(1,:)))