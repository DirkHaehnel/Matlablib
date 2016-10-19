if 0
    focpos = -4.9:0.1:0;
    for kk=1:21
        NA = 1.2 + 0.01*(kk-1);
        for j=1:length(focpos)
            [intx, inty, rho, phi] = ConfocalReflex([0 3], NA, 100, [1.52 1], [], 0.488, 5e3/1.3/1.8e3, focpos(j));
            col = cos(0*phi);
            subplot(5,10,j); 
            pcolor(cos(phi)*rho,sin(phi)*rho,intx+inty); 
            axis off; 
            axis image; 
            title(mnum2str(focpos(j),1,1),'fontname','times','fontsize',10);
            shading interp; 
            %        colormap gray
            drawnow
        end; 
        eval(['print -dpng -r300 oil' mint2str(100*NA,3) 'a']);
    end
    focpos = 0:0.1:4.9;
    for kk=1:21
        NA = 1.2 + 0.01*(kk-1);
        for j=1:length(focpos)
            [intx, inty, rho, phi] = ConfocalReflex([0 3], NA, 100, [1.52 1], [], 0.488, 5e3/1.3/1.8e3, focpos(j));
            col = cos(0*phi);
            subplot(5,10,j); 
            pcolor(cos(phi)*rho,sin(phi)*rho,intx+inty); 
            axis off; 
            axis image; 
            title(mnum2str(focpos(j),1,1),'fontname','times','fontsize',10);
            shading interp; 
            %        colormap gray
            drawnow
        end; 
        eval(['print -dpng -r300 oil' mint2str(100*NA,3) 'b']);
    end
end

if 0
    clear mask
    a = 199;
    b = 189;
    
    NAv = 1.21:0.01:1.3;
    fv = -0.35:0.05:-0.15;
    
    nn = 150;
    load D:\Joerg\Doc\Patra\ConfocalReflex\data.mat
    im = x(a-nn-10:a+nn+10,b-nn-10:b+nn+10,:);
    [x,y] = meshgrid(-nn:nn,-nn:nn);
    phi = angle(x+i*y);
    col = cos(0*phi);
    r = sqrt(x.^2+y.^2);
    bck = disk(nn);     
    pv = -5:0.5:5;
    for joerg = 1:length(pv)
        cnt = 1;
        for j=1:length(fv)
            for k=1:length(NAv)
                [fxx0, fxx2, byx0, byx2, rho] = ConfocalReflex([0 sqrt(2)*6.45*nn/100], NAv(k), 100, [1.52, 1.], [], 0.488, 5e3/1.3/1.8e3, pv(joerg)+fv(j));
                rho = rho/6.45;
                mask(:,:,cnt) = real((col.*interp1(rho,fxx0,r,'cubic')+cos(2*phi).*interp1(rho,fxx2,r,'cubic')).*...
                    conj(col.*interp1(rho,byx0,r,'cubic')+cos(2*phi).*interp1(rho,byx2,r,'cubic')))+...
                    real((sin(2*phi).*interp1(rho,fxx2,r,'cubic')).*conj(sin(2*phi).*interp1(rho,byx2,r,'cubic')));
                cnt = cnt+1
            end
        end
        for k=1:size(mask,3) mask(:,:,k) = mask(:,:,k)'.*bck; mask(:,:,k) = mask(:,:,k)/sum(sum(mask(:,:,k))); end
        eval(['save mask' mint2str(joerg,2) ' mask NAv fv'])
        [err, bim, cim, sim, xc{joerg}, yc{joerg}, bc{joerg}, cc{joerg}, sc{joerg}] = FindPattern(im(:,:,pos==pv(joerg)),mask,bck);
    end
    
    save temp
end

if 0
    clear mask
    NAv = 1.21:0.01:1.3;
    pv = -5:0.1:0;
    nn = 200;
    [x,y] = meshgrid(-nn:nn,-nn:nn);
    phi = angle(x+i*y);
    col = cos(0*phi);
    r = sqrt(x.^2+y.^2);
    bck = disk(nn);     
    for joerg = 1:length(pv)
        cnt = 1;
        for k=1:length(NAv)
            [fxx0, fxx2, byx0, byx2, rho] = ConfocalReflex([0 sqrt(2)*6.45*nn/100], NAv(k), 100, [1.52, 1.], [], 0.488, 5e3/1.3/1.8e3, pv(joerg));
            rho = rho/6.45;
            mask(:,:,cnt) = real((col.*interp1(rho,fxx0,r,'cubic')+cos(2*phi).*interp1(rho,fxx2,r,'cubic')).*...
                conj(col.*interp1(rho,byx0,r,'cubic')+cos(2*phi).*interp1(rho,byx2,r,'cubic')))+...
                real((sin(2*phi).*interp1(rho,fxx2,r,'cubic')).*conj(sin(2*phi).*interp1(rho,byx2,r,'cubic')));
            cnt = cnt+1
        end
        eval(['save mask' mint2str(10*pv(joerg),2) 'z mask NAv pv'])
    end
end

if 0
    load mask00z
    nn = (size(mask,1)-1)/2;
    [a,b]=meshgrid(-nn:nn,-nn:nn);
    r = sqrt(a.^2+b.^2);
    bin = 0.5:sqrt(2)*nn;
    for j=1:length(pv)
        eval(['load mask' mint2str(10*pv(j),2) 'z'])
        for k=1:length(NAv)
            tmp = mconv2(mask(:,:,k),disk(10));
            [ind ind] = max(tmp(:));
            res(j,k)=r(ind);
            %tmp = mask(:,:,k);
            %tmp = tmp>mconv2(tmp,disk(100));
            %h(:,j,k) = mhist(r(:),bin,tmp(:));
            %ind = bin(h(:,j,k)./bin'>2); 
            %res(j,k) = ind(end);
            %[ind ind]=min(gradient(h(:,j,k)./bin'));
            %res(j,k)=bin(ind);
        end
        j
    end
end

if 0
    load D:\Joerg\Doc\Patra\ConfocalReflex\data.mat 
    pv = -5:0.5:6;
    [a,b]=meshgrid(1:size(x,2),1:size(x,1));
    bin=0.5:1:300;
    for j=1:size(x,3)
        tmp=x(:,:,j)>mconv2(x(:,:,j),disk(100)); 
        ac(j)=sum(sum(a.*tmp))/sum(tmp(:)); 
        bc(j)=sum(sum(b.*tmp))/sum(tmp(:));    
        r=sqrt((a-ac(j)).^2+(b-bc(j)).^2);
        hh(:,j)=mhist(r(:),bin,tmp(:));
        ind=bin(hh(:,j)./bin'>4); 
        res(j)=ind(end);
    end
end