close all
clear all

prntbld = 0;

%name = 'D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2006-03-28\Water60xPH54\water60x.mat'
%name = 'D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2006-04-04\PSF_Water60x_180_180_540.mat'
%name = 'D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2006-03-01\PSFexp.mat' 
%name = 'D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2006-03-02\PSFexp.mat' 
%name = 'D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2006-03-02A\PSFexp.mat' 
name = 'D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2006-03-02B\PSFexp.mat' 
%name = 'D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2006-03-13\z_cc172_cs172_ph565_200um.mat' 
%name = 'D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2006-04-22\PSF.mat'
%name = 'D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2006-03-28\Oil\oil.mat' 
%name = 'D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_Iris_2006-08-16\PSFData.mat' 
%name = 'D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_Iris_2006-08-16\ZScanD.mat'; 
%name = 'D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_Iris_2006-08-17\ZScanA.mat'; 

eval(['load ' name]);

if 0 exa = ex; end
if 0 exa=exb; s = 'b'; else s = 'a'; end

if 1
    lamex = 640/1.33;
    lamem = 670/1.33;
else
    lamex = 640/1.51;
    lamem = 670/1.51;
end
av = 100/60*1e3;

if 1
    exa.wx1 = exa.wx1(2:end-1);
    exa.wy1 = exa.wy1(2:end-1);    
    exa.wx2 = exa.wx2(2:end-1);    
    exa.wy2 = exa.wy2(2:end-1);
	exa.z = exa.z(2:end-1);
	exa.amp1 = exa.amp1(2:end-1);
	exa.amp2 = exa.amp2(2:end-1);    
    exa.mx1 = exa.mx1(2:end-1);
    exa.my1 = exa.my1(2:end-1);    
    exa.mx2 = exa.mx2(2:end-1);    
    exa.my2 = exa.my2(2:end-1);
end
    
if 1
    facv = 1:0.01:1.15;
    z = 1e3*(exa.z/10-mean(exa.z/10));
    wx1 = exa.wx1; wx2 = exa.wx2;
    wy1 = exa.wy1; wy2 = exa.wy2;
    ind = abs(z)<=2000;
    for j=1:length(facv)
        fac = facv(j);
        l1 = simplex('LaserBeamFit',[200 100],[],[],[],[],lamex,z(ind),sqrt((wx1(ind).^2+wy1(ind).^2)/2)/fac,1);
        l2 = simplex('LaserBeamFit',l1,[],[],[],[],lamex,z(ind),sqrt((wx2(ind).^2+wy2(ind).^2)/2)/fac,1);
        err1(j) = LaserBeamFit(l1,lamex,z(ind),sqrt((wx1(ind).^2+wy1(ind).^2)/2)/fac,1);
        err2(j) = LaserBeamFit(l2,lamex,z(ind),sqrt((wx2(ind).^2+wy2(ind).^2)/2)/fac,1);
    end
    c1 = polyfit(facv,err1,2);
    c2 = polyfit(facv,err2,2);    
    [-c1(2)/c1(1)/2 -c2(2)/c2(1)/2]
    fac = mean([-c1(2)/c1(1)/2 -c2(2)/c2(1)/2])
else
    %fac = 1.0964
    fac = 1.019;
end
if 1    
    z = 1e3*(exa.z/10-mean(exa.z/10));
    wx1 = exa.wx1/fac; wx2 = exa.wx2/fac;
    wy1 = exa.wy1/fac; wy2 = exa.wy2/fac;
    %     wx1 = exa.wx1-50; wx2 = exa.wx2-50;
    %     wy1 = exa.wy1-50; wy2 = exa.wy2-50;
    amp1 = exa.amp1; amp2 = exa.amp2;
    dd = sqrt((exa.mx1-exa.mx2).^2+(exa.my1-exa.my2).^2)/fac;
end

if 1
    ind = abs(z)<=2000;
    close; l1=simplex('LaserBeamFit',[200 100],[],[],[],[],lamex,z(ind),sqrt((wx1(ind).^2+wy1(ind).^2)/2),1)
    close; l2=simplex('LaserBeamFit',l1,[],[],[],[],lamex,z(ind),sqrt((wx2(ind).^2+wy2(ind).^2)/2),1)
    [err,cl1] = LaserBeamFit(l1,lamex,z(ind),sqrt((wx1(ind).^2+wy1(ind).^2)/2),1);
    [err,cl2] = LaserBeamFit(l2,lamex,z(ind),sqrt((wx2(ind).^2+wy2(ind).^2)/2),1);
else
    close; l1=simplex('LaserBeamFit',[200 100],[],[],[],[],lamex,z,sqrt((wx1.^2+wy1.^2)/2))
    close; l2=simplex('LaserBeamFit',l1,[],[],[],[],lamex,z,sqrt((wx2.^2+wy2.^2)/2))
    [err,cl1] = LaserBeamFit(l1,lamex,z,sqrt((wx1.^2+wy1.^2)/2));
    [err,cl2] = LaserBeamFit(l2,lamex,z,sqrt((wx2.^2+wy2.^2)/2));
end
zi = z(1):z(end);
x1 = cl1*LaserBeamFit(l1,lamex,zi);
x2 = cl2*LaserBeamFit(l2,lamex,zi);

close; a1=simplex('D1CEF',[100 10],[],[],[],[],av,lamem,z,amp1,1)
close; a2=simplex('D1CEF',a1./[1.2 1]',[],[],[],[],av,lamem,z,amp2,1)
[err,c1]=D1CEF(a1,av,lamem,z,amp1,1);
[err,c2]=D1CEF(a2,av,lamem,z,amp2,1);
y1 = c1*D1CEF(a1,av,lamem,zi);
y2 = c2*D1CEF(a2,av,lamem,zi);
plot(z,amp1,'o',z,amp2,'o',zi,y1,'r',zi,y2,'b')
xlabel('\itz\rm-position [nm]')
ylabel('detection efficiency amplitude [a.u.]')
legend({'focus 1','focus 2'})
text(0.4,0.2,{['\ita\rm_1 = ' int2str(a1(1)) ' nm'],'',['\ita\rm_2 = ' int2str(a2(1)) ' nm']},'units','normal')

if prntbld
    eval(['print -dpng -r300 ' name(1:end-4) s '-Pinhole'])
end
    
figure

plot(z,sqrt((wx1.^2+wy1.^2)/2),'o',z,sqrt((wx2.^2+wy2.^2)/2),'o',zi,x1,'r',zi,x2,'b')
xlabel('\itz\rm-position [nm]')
ylabel('beam waist [nm]')
legend({'focus 1','focus 2'},3)
text(0.4,0.8,{['\itw\rm_1 = ' int2str(l1(1)) ' nm'],'',['\itw\rm_2 = ' int2str(l2(1)) ' nm']},'units','normal')

if prntbld
    eval(['print -dpng -r300 ' name(1:end-4) s '-Waist'])
end

vv = [GaussDetectionVolume(l1(1), a1(1), av, [lamex lamem]) GaussDetectionVolume(l2(1), a2(1), av, [lamex lamem])]*1e-9

figure
plot(z,dd); grid
xlabel('\itz\rm-position [nm]')
ylabel('distance [nm]')

if prntbld
    eval(['print -dpng -r300 ' name(1:end-4) s '-Distance'])
end

return

for j=1:20 
    p1(:,j) = simplex('LaserBeamFit',[200 100],[],[],[],[],lamex,z(ind),sqrt(((wx1(ind)-(j-1)*5).^2+(wy1(ind)-(j-1)*5).^2)/2),1,1);
    err1(j) = LaserBeamFit(p1(:,j),lamex,z(ind),sqrt(((wx1(ind)-(j-1)*5).^2+(wy1(ind)-(j-1)*5).^2)/2),1);
    p2(:,j) = simplex('LaserBeamFit',[200 100],[],[],[],[],lamex,z(ind),sqrt(((wx2(ind)-(j-1)*5).^2+(wy2(ind)-(j-1)*5).^2)/2),1,1);
    err2(j) = LaserBeamFit(p2(:,j),lamex,z(ind),sqrt(((wx2(ind)-(j-1)*5).^2+(wy2(ind)-(j-1)*5).^2)/2),1);
end

