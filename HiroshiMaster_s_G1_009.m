% HiroshiMaster

clear all

fname1 = 'D:\Daten\Johan\s_G1_008_1-500.TIF';
fname2 = 'D:\Daten\Johan\s_G1_008_501-1000.TIF'; 
outfile = 'D:\Daten\Johan\s_G1_008';

if 0
    filename = fname1;
    imnum = length(imfinfo(filename));
    % localize regions of interest
    im = double(imread(filename,1));
    for j=2:imnum
        im = im + double(imread(filename,j));
    end
    mim(im);
    pos = ginput;
    pos = round(pos);
end

pos = [...
    79    56
   206    30
   251    61
   300    45
   399    85
   304    97
    85   102
    68   137
   346   128
   318   141
   141   171
   165   194
   198   175
   223   209
   173   252
   142   268
    82   235
   165   299
   352   417
   258   396
   240   365
   437   269
   487   342
   420   329
   154   482
    54   189
    71   260
   131   328
    65   325
    60   386
   111   403
    91   451
   234   261
   249   202
   239   179
   281   292
   300   427
   306   242
   338   226
   365   267
   368   225
   373   191
   401   221
   438   120
   395   159
   436   233
   419   362
   406   307
   465   295
   404   468
   137   378
   147   382
   430   398
   418   481];

nn = 20;

if 0
    j=1;
    [res, z, zz, a, b, model] = Hiroshi(fname1,0,[pos(j,2)-nn pos(j,1)-nn; pos(j,2)+nn pos(j,1)+nn]);
    save HiroshiModel model
    eval(['save ' outfile 'M' mint2str(j,3) ' res z zz a b']);
    [res2, z2, zz2] = Hiroshi(fname2,0,[pos(j,2)-nn pos(j,1)-nn; pos(j,2)+nn pos(j,1)+nn],model);
    res = [res; res2];
    res(:,1) = res(:,1)+500*cumsum(diff([0; res(:,1)])<0);
    z = cat(3,z,z2);
    zz = cat(3,zz,zz2);
    eval(['save ' outfile 'M' mint2str(j,3) ' res z zz -append']);
    cnt = 1; bld = 1;
    for jj=1:size(z,3)
        subplot(4,5,cnt+5+5*floor((cnt-1)/5));
        mim(zz(:,:,jj));
        if ismember(jj,res(:,1))
            title({['\theta = ' mint2str(res(res(:,1)==jj,4)) '°']; ...
                ['\phi = ' mint2str(res(res(:,1)==jj,5)) '°']},'fontsize',12,'fontname','times');
        end
        subplot(4,5,cnt+5*floor((cnt-1)/5));
        mim(z(:,:,jj));
        title(num2str(jj),'fontsize',12);
        cnt = cnt+1;
        if cnt==11 | jj==ceil(size(z,3))
            eval(['print -dpng ' outfile 'M' mint2str(j,3) '-' mint2str(bld,3)]);
            bld = bld+1;
            cnt = 1;
        end
    end
else
    load HiroshiModel
end

for j=16:size(pos,1)
    [res, z, zz, a, b] = Hiroshi(fname1,0,[pos(j,2)-nn pos(j,1)-nn; pos(j,2)+nn pos(j,1)+nn],model);
    eval(['save ' outfile 'M' mint2str(j,3) ' res z zz a b']);
    [res2, z2, zz2] = Hiroshi(fname2,0,[pos(j,2)-nn pos(j,1)-nn; pos(j,2)+nn pos(j,1)+nn],model);
    res = [res; res2];
    res(:,1) = res(:,1)+500*cumsum(diff([0; res(:,1)])<0);
    z = cat(3,z,z2);
    zz = cat(3,zz,zz2);
    eval(['save ' outfile 'M' mint2str(j,3) ' res z zz -append']);
    cnt = 1; bld = 1;
    for jj=1:size(z,3)
        subplot(4,5,cnt+5+5*floor((cnt-1)/5));
        mim(zz(:,:,jj));
        if ismember(jj,res(:,1))
            title({['\theta = ' mint2str(res(res(:,1)==jj,4)) '°']; ...
                ['\phi = ' mint2str(res(res(:,1)==jj,5)) '°']},'fontsize',12,'fontname','times');
        end
        subplot(4,5,cnt+5*floor((cnt-1)/5));
        mim(z(:,:,jj));
        title(num2str(jj),'fontsize',12);
        cnt = cnt+1;
        if cnt==11 | jj==ceil(size(z,3))
            eval(['print -dpng ' outfile 'M' mint2str(j,3) '-' mint2str(bld,3)]);
            bld = bld+1;
            cnt = 1;
        end
    end
end

return

mim(im); hold on; h=plotsym(pos(:,1),pos(:,2),int2str((1:size(pos,1))')); hold off; for j=1:length(h) set(h(j),'fontsize',10,'color','y'); end

clear hh
bin = 0:0.5:50;
for k=1:size(pp,2)
    %eval(['load D:\Daten\Johan\s_G1_008M' mint2str(k,3)])
    eval(['load E:\s_G1_008M' mint2str(k,3)])
    theta=(90-abs(90-res(:,4)))/180*pi; phi=unwrap(res(:,5)/180*2*pi)/2;
    t = 1:size(res,1);
    t = t(~(res([1:end-1],1)+1==res([2:end],1)));
    t = [0 t size(res,1)];
    for j=1:length(t)-1
        pp(j,k).theta = theta(t(j)+1:t(j+1));
        pp(j,k).phi = phi(t(j)+1:t(j+1));
    end
    j = 1;
    tmp = [sin(pp(j,k).theta).*cos(pp(j,k).phi), sin(pp(j,k).theta).*sin(pp(j,k).phi), cos(pp(j,k).theta)];
    hh(:,k) = mhist(acos(sum(tmp(1:end-1,:).*tmp(2:end,:),2))/pi*180,bin);
    for j=2:size(pp,1)
        tmp = [sin(pp(j,k).theta).*cos(pp(j,k).phi), sin(pp(j,k).theta).*sin(pp(j,k).phi), cos(pp(j,k).theta)];
        hh(:,k) = hh(:,k) + mhist(acos(sum(tmp(1:end-1,:).*tmp(2:end,:),2))/pi*180,bin);
    end
end
plot(bin,hh./(ones(size(hh,1),1)*sum(hh)))


% prints time traces of the angles theta and phi
for k=1:size(pp,2)
    cnt = 0;
    j = 1;
    plot(cnt+1:cnt+length(pp(j,k).theta),pp(j,k).theta/pi*180,cnt+1:cnt+length(pp(j,k).theta),pp(j,k).phi/pi*180);
    cnt = cnt+length(pp(j,k).theta);
    hold on
    for j=2:size(pp,1)
        if ~isempty(pp(j,k).theta)
            plot(cnt+1:cnt+length(pp(j,k).theta),pp(j,k).theta/pi*180,cnt+1:cnt+length(pp(j,k).theta),pp(j,k).phi/pi*180);
            cnt = cnt+length(pp(j,k).theta);
        end
    end
    hold off
    xlabel('# of frame'); ylabel('angle [°]'); legend({'\theta','\phi'},5)
    eval(['print -dpng Trace' mint2str(k,3)])
end


% prints spherical projections of orientation traces 
for k=1:size(pp,2)
    sphere(100); shading interp;
    h=line([0 0],[0 0],[0 1.]); set(h,'color','c','linewidth',1);
    h=line([0 1.],[0 0],[0 0]); set(h,'color','c','linewidth',1);
    h=line([0 0],[0 0],[1.0 1.2]); set(h,'color','b','linewidth',1);
    h=line([1.0 1.2],[0 0],[0 0]); set(h,'color','b','linewidth',1);
    alpha(0.2); camlight(120, 60); axis image; axis off;
    hold on
    j = 1;
    plot3(sin(pp(j,k).theta).*cos(pp(j,k).phi), sin(pp(j,k).theta).*sin(pp(j,k).phi), cos(pp(j,k).theta),'linewidth',1);
    for j=2:size(pp,1)
        if ~isempty(pp(j,k).theta)
            plot3(sin(pp(j,k).theta).*cos(pp(j,k).phi), sin(pp(j,k).theta).*sin(pp(j,k).phi), cos(pp(j,k).theta),'linewidth',1);
        end
    end
    hold off
    view([120 30])
    axis([-1.2 1.2 -1.2 1.2 -1.2 1.2])
    axis vis3d
    eval(['print -dpng SphereTrace' mint2str(k,3)])
end


% in- and out-of-plane diffusion
bin = -25.125:0.25:25;
for k=1:size(pp,2)
    j = 1;
    tmp = [sin(pp(j,k).theta).*cos(pp(j,k).phi), sin(pp(j,k).theta).*sin(pp(j,k).phi), cos(pp(j,k).theta)];
    cc = cross(tmp(1:end-2,:),tmp(2:end-1,:),2);
    cc = cc./(sqrt(sum(cc.*cc,2))*[1 1 1]);
    dd = tmp(3:end,:)-tmp(2:end-1,:);
    hx(:,k) = mhist(asin(sum(cc.*dd,2))/pi*180,bin);
    cc = cross(cc,tmp(2:end-1,:),2);
    cc = cc./(sqrt(sum(cc.*cc,2))*[1 1 1]);
    hp(:,k) = mhist(asin(sum(cc.*dd,2))/pi*180,bin);
    for j=2:size(pp,1)
        if ~isempty(pp(j,k).theta)
            tmp = [sin(pp(j,k).theta).*cos(pp(j,k).phi), sin(pp(j,k).theta).*sin(pp(j,k).phi), cos(pp(j,k).theta)];
            cc = cross(tmp(1:end-2,:),tmp(2:end-1,:),2);
            cc = cc./(sqrt(sum(cc.*cc,2))*[1 1 1]);
            dd = tmp(3:end,:)-tmp(2:end-1,:);
            hx(:,k) = hx(:,k) + mhist(asin(sum(cc.*dd,2))/pi*180,bin);
            cc = cross(cc,tmp(2:end-1,:),2);
            cc = cc./(sqrt(sum(cc.*cc,2))*[1 1 1]);
            hp(:,k) = hp(:,k) + mhist(asin(sum(cc.*dd,2))/pi*180,bin);
        end
    end
end
plot(bin,sum(hx,2)/sum(hx(:)),bin,sum(hp,2)/sum(hp(:)))
for j=1:size(hx,2) plot(bin,hx(:,j)/sum(hx(:,j)),bin,hp(:,j)/sum(hp(:,j))); drawnow; pause; end


% Lorenz distribution fitting
for j=1:size(hx,2) 
    gx(:,j)=simplex('Lorenz',[0 1],[-inf 0],[],[],[],bin,hx(:,j),0); 
    [ex(j), cx(:,j)] = Lorenz(gx(:,j),bin,hx(:,j),0); 
    xlabel('out-of-plane angle change [°]'); ylabel('frequency')
    eval(['print -dpng xLorenz' mint2str(j,3)]) 
    gp(:,j)=simplex('Lorenz',[0 1],[-inf 0],[],[],[],bin,hp(:,j),0); 
    [ep(j), cp(:,j)] = Lorenz(gp(:,j),bin,hp(:,j),0); 
    xlabel('in-plane angle change [°]'); ylabel('frequency')
    eval(['print -dpng pLorenz' mint2str(j,3)])
end


% polynomial fitting of Lorenz width
c=polyfit([-gx(2,:) gx(2,:)],[-gp(2,:) gp(2,:)],3)
plot(gx(2,:),gp(2,:),'o',0:0.1:4.5,polyval(c,0:0.1:4.5))
xlabel('out-of-plane diffusion');ylabel('in-plane diffusion')


% Gaussian fitting
for j=1:size(hx,2) 
    gx(:,j)=simplex('Gauss',[0 0 1 2],[-inf -inf 0 0],[],[],[],bin,hx(:,j),0); 
    [ex(j), cx(:,j)] = Gauss(gx(:,j),bin,hx(:,j),0); 
    gp(:,j)=simplex('Gauss',[0 0 1 2],[-inf -inf 0 0],[],[],[],bin,hp(:,j),0); 
    [ep(j), cp(:,j)] = Gauss(gp(:,j),bin,hp(:,j),0); 
end


% Autocorrelation
kv = []; for j=1:size(pp,2) if length(pp(1,j).theta)==1e3 kv = [kv j]; end; end
clear auto
for cnt=1:length(kv)
    k=kv(cnt);
    j = 1;
    tmp = [sin(pp(j,k).theta).*cos(pp(j,k).phi), sin(pp(j,k).theta).*sin(pp(j,k).phi), cos(pp(j,k).theta)];
    for m=1:200 
        auto(m,cnt) = mean(acos(sum(tmp(1:end-m,:).*tmp(m+1:end,:),2)))/(size(tmp,1)-m-1); 
    end
end
plot(1:200,auto)
