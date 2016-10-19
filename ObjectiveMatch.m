% ObjectiveMatch

i = sqrt(-1);
mag = 100;
lambda = 0.65;

clear mask
nn = 40;
n2 = 2*nn;
n3 = (2*nn+1)^2;
[x,y] = meshgrid(-nn:nn,-nn:nn);
p = angle(x+i*y);
r = sqrt(x.^2+y.^2);

if 1
    
%     load 'C:\Daten\Digambara\NA\NewObjective\QD650_h9_water10000times_1,3objOy_15mW.mat'
%     names = {'h9a', 'h9b', 'h9c', 'h9d'};
%     NAv = 1.23;
%     abv = -0.02;
%     dzv = 0.65;
%     mm = 5;

  
    load 'C:\Daten\Digambara\NA\NewObjective\QD650_z1_water10000times_1,3objOy_15mW.mat' 
    names = {'z1a', 'z1b', 'z1c'};
    
    NAv = 1.22;
    abv = -0.02;
    dzv = 0.2;
    mm = 9;
    
    for jf=1:8
        for jNA=1:length(NAv)
            focus = (jf-1)*0.5 - dzv(jNA);
            NA = NAv(jNA);
            alpha = [zeros(1,mm) 2*pi*abv(jNA)];
            
            [intx inty intz, rho, phi, fxx0, fxx2, fxz, byx0, byx2, byz] = ...
                SEPDipole([0 max(r(:))*6.45/mag], 0, NA, 1.52, 1.52, 1.0, [], 0, [], lambda, mag, focus, [], alpha);
            rho = rho/6.45;
            cnt = 1;
            for j=1:16
                al = 0; be = pi/16*(j-1);
                int = real((cos(al)*(interp1(rho,fxx0,r,'cubic')+cos(2*(p-be)).*interp1(rho,fxx2,r,'cubic'))+sin(al)*cos(p-be).*interp1(rho,fxz,r,'cubic')).*...
                    conj(cos(al)*(interp1(rho,byx0,r,'cubic')+cos(2*(p-be)).*interp1(rho,byx2,r,'cubic'))+sin(al)*cos(p-be).*interp1(rho,byz,r,'cubic')) + ...
                    (cos(al)*sin(2*(p-be)).*interp1(rho,fxx2,r,'cubic')+sin(al)*sin(p-be).*interp1(rho,fxz,r,'cubic')).*...
                    conj(cos(al)*sin(2*(p-be)).*interp1(rho,byx2,r,'cubic')+sin(al)*sin(p-be).*interp1(rho,byz,r,'cubic')));    
                mask(:,:,cnt) = int; 
                theta(cnt) = al; phi(cnt) = be;
                cnt = cnt+1;
            end
            
            for j=1:size(mask,3) 
                subplot(3,6,j); 
                mim(mask(:,:,j)); 
                line([0.5 size(mask,2)+0.5 size(mask,2)+0.5 0.5 0.5],[0.5 0.5 size(mask,1)+0.5 size(mask,1)+0.5 0.5],'color','k','linewidth',1)
                title(['\theta = ' mint2str(theta(j)/pi*180) '°, \phi = ' mint2str(phi(j)/pi*180) '°'],'fontsize',7,'fontname','times'); 
            end
            colormap(flipud(gray))
            bck = disk(nn); 
            for j=1:size(mask,3) mask(:,:,j) = mask(:,:,j).*bck; mask(:,:,j) = mask(:,:,j)/sum(sum(mask(:,:,j))); end
            
            for jn = 1:length(names)
                clear err bim cim sim xc yc z orient intens
                eval(['im = ' names{jn} '(:,:,jf);']);
                
                [err, bim, cim, sim, xc, yc] = FindPattern(im,mask,bck,0.7,1,'mconv2(cim,disk(7))>tsh*sqrt(err)',false);

                tmp = 0*im;
                [a,b] = size(tmp);
                
                for k=1:length(xc)
                    tmp(max(yc(k)-nn,1):min(yc(k)+nn,a),max(xc(k)-nn,1):min(xc(k)+nn,b)) = ...
                        tmp(max(yc(k)-nn,1):min(yc(k)+nn,a),max(xc(k)-nn,1):min(xc(k)+nn,b)) + ...
                        cim(yc(k),xc(k))*mask(max(nn-yc(k)+2,1):min(2*nn+1-yc(k)-nn+a,2*nn+1),...
                        max(nn-xc(k)+2,1):min(2*nn+1-xc(k)-nn+b,2*nn+1),sim(yc(k),xc(k)));
                end
                
                eval([names{jn} 'res(:,:,jf,jNA) = tmp;']);
            end
                
        end
    end
end

