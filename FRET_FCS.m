if 1 % FRET-FCS over glass surface in water 
    close all
    
    %read emission spec
    fin = fopen('c:\Joerg\Doc\Microscopy\Defocused\Diamond\NV Emission.dat','r');
    fgetl(fin);
    nvem = fscanf(fin,'%f %f\n');
    nvem = [nvem(1:2:end) nvem(2:2:end)];
    
    lamex = nvem(nvem(:,2)==max(nvem(:,2)),1);
    
    %read absorption spec
    fin = fopen('c:\Joerg\Doc\Microscopy\Defocused\Diamond\Abs-IR800.dat','r');
    fgetl(fin);
    irab = fscanf(fin,'%f %f\n');
    irab = [irab(1:2:end) irab(2:2:end)];
    
    lamem = irab(irab(:,2)==max(irab(:,2)),1);
    
    rad = 2*pi/lamex*0.01*exp(log(10/0.01)*(0:0.01:1));
    z = 2*pi/lamex*0.01*exp(log(10/0.01)*(0:0.01:1));
    zd = 0;
    n0 = 2.42;
    nd = 2.42;
    ni = 1.33;
    na = 1.33;
    n1 = 1.33;
    d0 = [];
    dd = 2*pi/lamex*5;
    di = 0;
    da = 0;
    d1 = [];

    ex = zeros(3,length(rad),length(z));
    ex0 = ex;
    ez = zeros(2,length(rad),length(z));
    ez0 = ez;
    
    for jz = 1:length(z)
        [ex(:,:,jz),ez(:,:,jz),ex0(:,:,jz),ez0(:,:,jz)] = DipoleField(zd,rad,z(jz),n0,nd,ni,na,n1,d0,dd,di,da,d1);
        jz
    end
    exctotal = squeeze(abs(ex0(1,:,:)).^2 + abs(ex0(2,:,:)).^2 + 0.5*abs(ex0(3,:,:)).^2 + abs(ez0(1,:,:)).^2 + abs(ez0(2,:,:)).^2);
    pcolor([-lamex/2/pi*rad(end:-1:1)'*ones(1,numel(z)); lamex/2/pi*rad'*ones(1,numel(z))],lamex/2/pi*ones(2*numel(rad),1)*z,log10([exctotal(end:-1:1,:); exctotal]/max(exctotal(:)))); shading interp; axis image; colorbar

    figure
    exctotal = squeeze(abs(ex(1,:,:)).^2 + abs(ex(2,:,:)).^2 + 0.5*abs(ex(3,:,:)).^2 + abs(ez(1,:,:)).^2 + abs(ez(2,:,:)).^2);
    pcolor([-lamex/2/pi*rad(end:-1:1)'*ones(1,numel(z)); lamex/2/pi*rad'*ones(1,numel(z))],lamex/2/pi*ones(2*numel(rad),1)*z,log10([exctotal(end:-1:1,:); exctotal]/max(exctotal(:)))); shading interp; axis image; colorbar

    figure 
    FocusImage3D(rad'*ones(1,numel(z)),ones(numel(rad),1)*z,cat(3,squeeze(ex(1,:,:)),zeros(size(ex,2),size(ex,3)),squeeze(ex(2,:,:)),zeros(size(ex,2),size(ex,3)),zeros(size(ex,2),size(ex,3))));
    
    %[modres, autotime] = FCS(lamex*rad/2e3/pi,lamex*z/2e3/pi,exctotal);
    
end

