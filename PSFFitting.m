%%% ERROR ANALYSIS %%%

%load D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2006-02-20\PSFexp.mat
%exa.z = al*exa.z; exb.z = exb.z/al-be;

%load D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2006-02-28\PSFexp.mat
%al = al1; be = be1; exa.z = al*exa.z*1e2 - 31e3; exb.z = (exb.z/al-be)*1e2 - 31e3;
%al = al2; be = be2; exa = exe; exb = exf; exa.z = al*exa.z*1e2 - 29.5e3; exb.z = (exb.z/al-be)*1e2 - 29.5e3;

% load D:\Doc\Fcs\2Focus\PSF\PSF_2006-03-01\PSFexp.mat
% exa.z = exa.z*1e2 - 31e3; exb.z = exb.z*1e2 - 31e3;

%load D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2006-03-02\PSFexp.mat
% al = 1.035; be=-(al-1/al)*max(exa.z); plot(al*exa.z,exa.wx1,exb.z/al-be,exb.wx1)
% exa.z = al*exa.z*1e2 - 31e3; exb.z = (exb.z/al-be)*1e2 - 31e3;

% load D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2006-03-02B\PSFexp.mat
% exa.z = exa.z*1e2 - 32e3; exb.z = exb.z*1e2 - 32e3;

load D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2006-03-13\z_cc172_cs172_ph565_200umJoerg.mat 
al = 1.012; be=-(al-1/al)*max(exa.z); plot(al*exa.z,exa.wx1,exb.z/al-be,exb.wx1)
% exa.z = al*exa.z*1e2 - 21e3; exb.z = (exb.z/al-be)*1e2 - 21e3;
exa.z = al*exa.z; exb.z = (exb.z/al-be);

% [pos, pos] = min(exa.wx1);
% pos = pos-1;
% tmp = [exa.wx1(1:pos) exa.wy1(pos+1:end)];
% exa.wy1 = [exa.wy1(1:pos) exa.wx1(pos+1:end)];
% exa.wx1 = tmp;
% [pos, pos] = min(exa.wx2);
% pos = pos - 1;
% tmp = [exa.wx2(1:pos) exa.wy2(pos+1:end)];
% exa.wy2 = [exa.wy2(1:pos) exa.wx2(pos+1:end)];
% exa.wx2 = tmp;
% [pos, pos] = min(exb.wx1);
% pos = pos - 1;
% tmp = [exb.wx1(1:pos) exb.wy1(pos+1:end)];
% exb.wy1 = [exb.wy1(1:pos) exb.wx1(pos+1:end)];
% exb.wx1 = tmp;
% [pos, pos] = min(exb.wx2);
% pos = pos - 1;
% tmp = [exb.wx2(1:pos) exb.wy2(pos+1:end)];
% exb.wy2 = [exb.wy2(1:pos) exb.wx2(pos+1:end)];
% exb.wx2 = tmp;

load PSFModel_NA114
z = -z;

lim = 2e3;
inda = exa.z>=-lim & exa.z<=lim;
indb = exb.z>=-lim & exb.z<=lim;
shiftv = (-1e3:10:1e3);
[vshift,vrad,vcov,vabe] = ndgrid(1:length(shiftv),1:length(radv),1:length(covv),1:length(abev));
clear err
for jrad = 1:length(radv)
    for jcov = 1:length(covv)
        for jabe = 1:length(abev)
            for jshift=1:length(shiftv)
                shift = shiftv(jshift);
                ind = imag(wx1(:,jrad,jcov,jabe))==0;
                wxa = interp1(shift+1e3*z(1,ind),1e3*wx1(ind,jrad,jcov,jabe),exa.z(inda));
                wxb = interp1(shift+1e3*z(1,ind),1e3*wx1(ind,jrad,jcov,jabe),exb.z(indb));
                err(jshift,jrad,jcov,jabe) = sum((exa.wx1(inda)-wxa).^2./wxa)+sum((exb.wx1(indb)-wxb).^2./wxb);
            end
        end
    end
end
jshiftx1 = vshift(err==min(err(:)));
jradx1 = vrad(err==min(err(:)));
jcovx1 = vcov(err==min(err(:)));
jabex1 = vabe(err==min(err(:)));
for jrad = 1:length(radv)
    for jcov = 1:length(covv)
        for jabe = 1:length(abev)
            for jshift=1:length(shiftv)
                shift = shiftv(jshift);
                ind = imag(wy1(:,jrad,jcov,jabe))==0;
                wya = interp1(shift+1e3*z(1,ind),1e3*wy1(ind,jrad,jcov,jabe),exa.z(inda));
                wyb = interp1(shift+1e3*z(1,ind),1e3*wy1(ind,jrad,jcov,jabe),exb.z(indb));
                err(jshift,jrad,jcov,jabe) = sum((exa.wy1(inda)-wya).^2./wya)+sum((exb.wy1(indb)-wyb).^2./wyb);
            end
        end
    end
end
jshifty1 = vshift(err==min(err(:)));
jrady1 = vrad(err==min(err(:)));
jcovy1 = vcov(err==min(err(:)));
jabey1 = vabe(err==min(err(:)));

plot(exa.z,exa.wx1,'o',exa.z,exa.wy1,'o',exb.z,exb.wx1,'x',exb.z,exb.wy1,'x',shiftv(jshiftx1)+1e3*z(1,:),1e3*wx1(:,jradx1,jcovx1,jabex1),shiftv(jshifty1)+1e3*z(1,:),1e3*wy1(:,jrady1,jcovy1,jabey1))
xlabel('z-position [nm]'); ylabel('beam waists of 1st focus [nm]');


clear err
for jrad = 1:length(radv)
    for jcov = 1:length(covv)
        for jabe = 1:length(abev)
            for jshift=1:length(shiftv)
                shift = shiftv(jshift);
                ind = imag(wx2(:,jrad,jcov,jabe))==0;
                wxa = interp1(shift+1e3*z(1,ind),1e3*wx2(ind,jrad,jcov,jabe),exa.z(inda));
                wxb = interp1(shift+1e3*z(1,ind),1e3*wx2(ind,jrad,jcov,jabe),exb.z(indb));
                err(jshift,jrad,jcov,jabe) = sum((exa.wx2(inda)-wxa).^2./wxa)+sum((exb.wx2(indb)-wxb).^2./wxb);
            end
        end
    end
end
jshiftx2 = vshift(err==min(err(:)));
jradx2 = vrad(err==min(err(:)));
jcovx2 = vcov(err==min(err(:)));
jabex2 = vabe(err==min(err(:)));
for jrad = 1:length(radv)
    for jcov = 1:length(covv)
        for jabe = 1:length(abev)
            for jshift=1:length(shiftv)
                shift = shiftv(jshift);
                ind = imag(wy2(:,jrad,jcov,jabe))==0;
                wya = interp1(shift+1e3*z(1,ind),1e3*wy2(ind,jrad,jcov,jabe),exa.z(inda));
                wyb = interp1(shift+1e3*z(1,ind),1e3*wy2(ind,jrad,jcov,jabe),exb.z(indb));
                err(jshift,jrad,jcov,jabe) = sum((exa.wy2(inda)-wya).^2./wya)+sum((exb.wy2(indb)-wyb).^2./wyb);
            end
        end
    end
end
jshifty2 = vshift(err==min(err(:)));
jrady2 = vrad(err==min(err(:)));
jcovy2 = vcov(err==min(err(:)));
jabey2 = vabe(err==min(err(:)));

plot(exa.z,exa.wx2,'o',exa.z,exa.wy2,'o',exb.z,exb.wx2,'x',exb.z,exb.wy2,'x',shiftv(jshiftx2)+1e3*z(1,:),1e3*wx2(:,jradx2,jcovx2,jabex2),shiftv(jshifty2)+1e3*z(1,:),1e3*wy2(:,jrady2,jcovy2,jabey2))
xlabel('z-position [nm]'); ylabel('beam waists of 1st focus [nm]');


shift=[shiftv(jshiftx1) shiftv(jshifty1) shiftv(jshiftx2) shiftv(jshifty2)]
rad=[radv(jradx1) radv(jrady1) radv(jradx2) radv(jrady2)]
cov=[covv(jcovx1) covv(jcovy1) covv(jcovx2) covv(jcovy2)]

da=sqrt((exa.mx1-exa.mx2).^2+(exa.my1-exa.my2).^2);
db=sqrt((exb.mx1-exb.mx2).^2+(exb.my1-exb.my2).^2);
ind1 = imag(mx1(:,jradx1,jcovx1))==0 & imag(my1(:,jrady1,jcovy1))==0;
ind2 = imag(mx2(:,jradx2,jcovx2))==0 & imag(my2(:,jrady2,jcovy2))==0;
dda = sqrt((interp1(shiftv(jshiftx1)+1e3*z(1,ind1),1e3*mx1(:,jradx1,jcovx1),exa.z)-...
    interp1(shiftv(jshiftx2)+1e3*z(1,ind2),1e3*mx2(:,jradx2,jcovx2),exa.z)).^2 + ...
    (interp1(shiftv(jshifty1)+1e3*z(1,ind1),1e3*my1(:,jrady1,jcovy1),exa.z)-...
    interp1(shiftv(jshifty2)+1e3*z(1,ind2),1e3*my2(:,jrady2,jcovy2),exa.z)).^2);
ddb = sqrt((interp1(shiftv(jshiftx1)+1e3*z(1,ind1),1e3*mx1(:,jradx1,jcovx1),exb.z)-...
    interp1(shiftv(jshiftx2)+1e3*z(1,ind2),1e3*mx2(:,jradx2,jcovx2),exb.z)).^2 + ...
    (interp1(shiftv(jshifty1)+1e3*z(1,ind1),1e3*my1(:,jrady1,jcovy1),exb.z)-...
    interp1(shiftv(jshifty2)+1e3*z(1,ind2),1e3*my2(:,jrady2,jcovy2),exb.z)).^2);
plot(exa.z,da,'o-',exa.z,dda,exb.z,db,'x-',exb.z,ddb)
xlabel('z-position [nm]'); ylabel('focus distance [nm]');

plot(exa.z,exa.amp1/max(exa.amp1),'o',exa.z,exa.amp2/max(exa.amp2),'o',...
    exb.z,exb.amp1/max(exb.amp1),'x',exb.z,exb.amp2/max(exb.amp2),'x',...
    shiftv(jshiftx1)+1e3*z(1,:),amp1(:,jradx1,jcovx1)/max(amp1(:,jradx1,jcovx1)),...
    shiftv(jshiftx2)+1e3*z(1,:),amp2(:,jradx2,jcovx2)/max(amp2(:,jradx2,jcovx2)),...
    shiftv(jshifty1)+1e3*z(1,:),amp1(:,jrady1,jcovy1)/max(amp1(:,jrady1,jcovy1)),...
    shiftv(jshifty2)+1e3*z(1,:),amp2(:,jrady2,jcovy2)/max(amp2(:,jrady2,jcovy2)))
xlabel('z-position [nm]'); ylabel('normalized amplitude');

return % Amplitude fitting

lim = 2e3;
inda = exa.z>=-lim & exa.z<=lim;
indb = exb.z>=-lim & exb.z<=lim;
shiftv = (-1e3:10:1e3);
[vshift,vrad,vcov] = ndgrid(1:length(shiftv),1:length(radv),1:length(covv));
clear err
for jrad = 1:length(radv)
    for jcov = 1:length(covv)
        for jshift=1:length(shiftv)
            shift = shiftv(jshift);
            ind = imag(wx1(:,jrad,jcov,jabe))==0;
            a = interp1(shift+1e3*z(1,ind),amp1(ind,jrad,jcov,jabe),exa.z(inda))';
            b = [a>0 a]*([a>0 a]\exb.amp1(inda)');
            a = [a>0 a]*([a>0 a]\exa.amp1(inda)');
            err(jshift,jrad,jcov,jabe) = sum((exa.amp1(inda)'-a).^2./a) + sum((exb.amp1(inda)'-b).^2./b);
        end
    end
end
jshift1 = vshift(err==min(err(:)));
jrad1 = vrad(err==min(err(:)));
jcov1 = vcov(err==min(err(:)));
for jrad = 1:length(radv)
    for jcov = 1:length(covv)
        for jshift=1:length(shiftv)
            shift = shiftv(jshift);
            ind = imag(wy1(:,jrad,jcov,jabe))==0;
            a = interp1(shift+1e3*z(1,ind),amp2(ind,jrad,jcov,jabe),exa.z(indb))';
            b = [a>0 a]*([a>0 a]\exb.amp2(indb)');
            a = [a>0 a]*([a>0 a]\exa.amp2(indb)');
            err(jshift,jrad,jcov,jabe) = sum((exa.amp2(indb)'-a).^2./a) + sum((exb.amp2(indb)'-b).^2./b);
        end
    end
end
jshift2 = vshift(err==min(err(:)));
jrad2 = vrad(err==min(err(:)));
jcov2 = vcov(err==min(err(:)));

a = interp1(shiftv(jshift1)+1e3*z(1,ind),amp1(ind,jrad1,jcov1),exa.z(inda))';
a = [ones(size(z,2),1) amp1(:,jrad1,jcov1)]*([a>0 a]\exa.amp1(inda)');
b = interp1(shiftv(jshift2)+1e3*z(1,ind),amp2(ind,jrad2,jcov2),exa.z(indb))';
b = [ones(size(z,2),1) amp2(:,jrad1,jcov1)]*([b>0 b]\exb.amp2(indb)');

plot(exa.z,exa.amp1,'o',exa.z,exa.amp2,'x',shiftv(jshift1)+1e3*z(1,:),a,shiftv(jshift2)+1e3*z(1,:),b)
xlabel('z-position [nm]'); ylabel('amplitude');


return

load D:\Joerg\Doc\Fcs\2Focus\PSF\PSF_2006-03-28\Bead1.mat

clear ex mx1 mx2 my1 my2 wx1 wx2 wy1 wy2 amp1 amp2 d 
im1 = exb.im1;
im2 = exb.im2;
z = exb.z;

pix = 200e3/2^12;
for j=1:size(im1,3)
	[mx1(j),my1(j),wx1(j),wy1(j),amp1(j)] = Gauss2D([],[],im1(:,:,j),pi/4);
    [mx2(j),my2(j),wx2(j),wy2(j),amp2(j)] = Gauss2D([],[],im2(:,:,j),pi/4);
    d(j) = sqrt((mx1(j)-mx2(j))^2+(my1(j)-my2(j))^2)*pix;
end
ex.wx1 = pix*wx1;
ex.wx2 = pix*wx2;
ex.wy1 = pix*wy1;
ex.wy2 = pix*wy2;
ex.mx1 = pix*mx1;
ex.mx2 = pix*mx2;
ex.my1 = pix*my1;
ex.my2 = pix*my2;
ex.amp1 = amp1;
ex.amp2 = amp2;
ex.im1 = im1;
ex.im2 = im2;
ex.z = z;
ex.d = d;
