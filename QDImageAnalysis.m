% QD Image Analysis

if 0
    load D:\Joerg\Doc\Patra\SauerDots\picks.mat    
    load QD2d1d
    
    j=1; h = QDControl('Image', imm(:,:,j), 'theta', res2d1d(j,1), 'omega', res2d1d(j,2), 'phi', res2d1d(j,3), 'kappa', res2d1d(j,4), 'ratio', res2d1d(j,5), 'focus', res2d1d(j,6))
    tmp = guidata(h); res3(j,:) = [tmp.theta tmp.omega tmp.phi tmp.kappa tmp.ratio tmp.focus tmp.mag tmp.na]; pic3(:,:,j) = tmp.int; save QDtriple res3 pic3 -append
    
    % load D:\Joerg\Doc\Patra\SauerDots\picks.mat    
    % load QDtriple
    % ind = [3 19 16 36 41 23 14 22 26 31];
    % for j=1:5 subplot(5,5,j); mim(imm(:,:,ind(j))); subplot(5,5,j+5); mim(pic3(:,:,ind(j))); subplot(5,5,j+15); mim(imm(:,:,ind(j+5))); subplot(5,5,j+20); mim(pic3(:,:,ind(j+5))); end
    % ord = [5 10 4 9 3 8 2 7 1 6];
    % for j=1:10 tab(j,:) = res3(ind(ord(j)),1:end-2); end
    % tab(:,3) = mod(360-tab(:,3)+90+180,360)-180
    % fout=fopen('Fig4Param.txt','w'); fprintf(fout,'%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\n',tab'); fclose(fout)
end

if 0
    NA = 1.2;
    mag = 110;
    lambda = 0.57;
    
    clear mask
    nn = 18;
    n2 = 2*nn;
    n3 = (2*nn+1)^2;
    [x,y] = meshgrid(-nn:nn,-nn:nn);
    p = angle(x+i*y);
    r = sqrt(x.^2+y.^2);
    
    cnt = 1;
    for jfoc=1:5
        focus = 1 + 0.05*(jfoc-1);        
        for jom=1:7
            om = pi/12*(jom-1);
            [tmp, rho, phi, intc1, ints1] = QDSEP(0,0,om,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus);
            [tmp, rho, phi, intc2, ints2] = QDSEP(0,pi/2,om,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus,pi/2);
            rho = rho/6.45;
            int = interp1(rho,real(intc1(1,:)+intc2(1,:)),r,'cubic');
            for jj=1:4 
                int = int + cos(jj*p).*interp1(rho,real(intc1(jj+1,:)),r,'cubic') + ...
                    sin(jj*p).*interp1(rho,real(ints1(jj,:)),r,'cubic'); 
                int = int + cos(jj*(p+pi/2)).*interp1(rho,real(intc2(jj+1,:)),r,'cubic') + ...
                    sin(jj*(p+pi/2)).*interp1(rho,real(ints2(jj,:)),r,'cubic');                 
            end
            [tmp, tmp, tmp, intc1, ints1] = QDSEP(-1,0,om+pi/2,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus);
            [tmp, tmp, tmp, intc2, ints2] = QDSEP(-1,pi/2,om+pi/2,0,2,0,NA,1.51,1,1,[],0,[],lambda,mag,focus,pi/2);
            tmp = interp1(rho,real(intc1(1,:)+intc2(1,:)),r,'cubic');
            for jj=1:4 
                tmp = tmp + cos(jj*p).*interp1(rho,real(intc1(jj+1,:)),r,'cubic') + ...
                    sin(jj*p).*interp1(rho,real(ints1(jj,:)),r,'cubic'); 
                tmp = tmp + cos(jj*(p+pi/2)).*interp1(rho,real(intc2(jj+1,:)),r,'cubic') + ...
                    sin(jj*(p+pi/2)).*interp1(rho,real(ints2(jj,:)),r,'cubic');                 
            end
            
            for jkappa=1:7
                kappa = (jkappa-1)/6;
                mask(:,:,cnt) = (1-kappa)*int + kappa*tmp;
                kappav(cnt) = kappa;
                omv(cnt) = om;
                focv(cnt) = focus;
                cnt = cnt + 1;
            end
        end
        save maskres mask kappav omv focv
    end
end

if 0
    load maskres
    orient landscape
    for jfoc = 2:5
        focus = 1 + 0.05*(jfoc-1);
        clf 
        cnt = 1;
        for jom=1:7
            om = pi/12*(jom-1);
            for jkappa=1:7
                kappa = (jkappa-1)/6;
                subplot(7,7,cnt)
                mim(mask(:,:,(jfoc-1)*7*7+(jom-1)*7+jkappa));
                if jom==1
                    title(mnum2str(kappa,1,2),'fontname','times','fontsize',10);
                end
                if jkappa==1
                    h = text(-10,21,int2str(om/pi*180),'fontname','times','fontsize',10);
                    set(h,'horizontalalignment','right');
                end
                drawnow
                cnt = cnt+1;
            end
        end
        eval(['print -dpng -r300 focus' mint2str(100*focus,3)]);        
    end
end

if 0
    load D:\Joerg\Doc\Patra\SauerDots\picks.mat
    load D:\Joerg\Doc\Patra\SauerDots\pickangles.mat
    load maskres
    n1 = size(mask,1);
    n2 = size(mask,2);
    n3 = size(mask,3);
    bck = ones(n1,n2);
    bck = bck/sqrt(sum(sum(bck.^2)));
    for j=1:n3 
        mask(:,:,j) = mask(:,:,j).*(bck>0);
        mask(:,:,j) = mask(:,:,j)./sqrt(sum(sum(mask(:,:,j).^2)));
    end
    nn1 = size(imm,1);
    nn2 = size(imm,2);
    nn3 = size(imm,3);
    
    for k=1:nn3
        err = inf*ones(nn1,nn2);
        cim = ones(nn1,nn2);
        bim = ones(nn1,nn2);
        sim = ones(nn1,nn2);
        
        im = imrotate(imm(:,:,k),phi(k),'bicubic','crop');
        for s = 1:n3
            im0 = mconv2(im,bck);
            im02 = mconv2(im.^2,bck/max(max(bck)));
            crs = sum(sum(bck.*mask(:,:,s)));
            crs = inv([1 crs; crs 1]);
            im1 = mconv2(im,squeeze(mask(:,:,s)));
            tmperr = im02 - crs(1,1)*im0.*im0 - 2*crs(1,2)*im0.*im1 - crs(2,2)*im1.*im1;
            tmpbim = crs(1,1)*im0 + crs(1,2)*im1;
            tmpcim = crs(1,2)*im0 + crs(2,2)*im1;
            cim(tmperr<err) = tmpcim(tmperr<err);
            bim(tmperr<err) = tmpbim(tmperr<err);
            sim(tmperr<err) = s;
            err(tmperr<err) = tmperr(tmperr<err);
            s
        end
        cres(:,:,k) = cim;
        bres(:,:,k) = bim;
        sres(:,:,k) = sim;
        eres(:,:,k) = err;
        save QDIAres imm cres bres sres eres
    end
end

break

for j=1:size(imm,3)
    ind = (eres(19:23,19:23,j)./sqrt(cres(19:23,19:23,j)))==min(min(eres(19:23,19:23,j)./sqrt(cres(19:23,19:23,j))));
    subplot(121); 
    tmp = rres(19:23,19:23,j);
    mim(imrotate(imm(:,:,j),10*(tmp(ind)-1),'bicubic','crop')); 
    subplot(122); 
    tmp = sres(19:23,19:23,j);
    mim(mask(:,:,tmp(ind)));
    pause
end

break 

for j=1:size(imm,3) 
    tst=sres(21-3:21+3,21-3:21+3,j); 
    tst=tst(eres(21-3:21+3,21-3:21+3,j)==min(min(eres(21-3:21+3,21-3:21+3,j))));
    subplot(121); mim(imrotate(imm(:,:,j),phi(j),'bicubic','crop')); subplot(122); tmp=[zeros(4,41); [zeros(33,4) mask(:,:,tst) zeros(33,4)]; zeros(4,41)]; mim(tmp);
    pause
end

break 

if 0
    load maskres
    orient landscape
    for jfoc = 1:5
        focus = 1 + 0.05*(jfoc-1);
        for panel = 1:4
            clf 
            cnt = 1;
            for jom=(rem(2,panel)*5+(1:11))
                om = pi/40*(jom-1);
                for jkappa=(mod(panel-1,2)*10+(1:11))
                    kappa = jkappa/20;
                    subplot(11,11,cnt)
                    mim(mask(:,:,(jom-1)*7*21+(jfoc-1)*21+jkappa));
                    if ((panel==1 | panel==2) & jom==1) | ((panel==3 | panel==4) & jom==11)
                        title(mnum2str(kappa,1,2),'fontname','times','fontsize',10);
                    end
                    if ((panel==1 | panel==3) & jkappa==1) | ((panel==2 | panel==4) & jkappa==11)
                        h = text(-10,21,int2str(om/pi*180),'fontname','times','fontsize',10);
                        set(h,'horizontalalignment','right');
                    end
                    drawnow
                    cnt = cnt+1;
                end
            end
            eval(['print -dpng -r300 focus' mint2str(100*focus,3) 'panel' mint2str(panel,1)]);
            print
        end
    end
end
