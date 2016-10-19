function [err, c, z, zz] = Fcs1D(p,delta,y,t,expflag,vflag,pic)
% function [err, c, z, zz] = Fcs1D(p,delta,y,t,expflag,vflag,pic)
% Is used for fitting of 2fFCS diffusion and flow data for Nanochannels. 
%
% Input Parameters:
%
% p(1)                      - Beam waist radius for focus 1
% p(2)                      - Beam waist radius for focus 2
% p(3:end-vflag-expflag)    - diffusion coefficient (um^2/s)
% p(numel(pd)+3:end-vflag)  - exponential components (triplet states)
% p(end-vflag+1:end)        - velocity value
% delta                     - calibrated distance between the two focii
% y                         - data points
% t                         - lag time
% expflag                   - Number of exponential components
% vflag                     - Number of velocity components
% pic                       - if 1, shows the fitting optimization with the
%                             data points.
%
% Usage: p = Simplex('Fcs1D',p,[],[],[],1000,delta,y,t,1,1,1)
% Refer to Dertinger et al., Optics Express, 16, 14353 (2008)
% or type on the link : http://www.joerg-enderlein.de/fileadmin/publications/oe08c.pdf

% (c) Narain Karedla (2014)



if numel(size(y))>2
    y=mean(y,3);
end
if nargin<6 || isempty(pic)
    pic=0;
end
if nargin<5 || isempty(vflag)
    vflag=0;
end

if nargin<4 || isempty(expflag)
    expflag=0;
end
% p=1/p;
w1       = p(1); % beam waist radius for focus 1
w2       = p(2); % beam waist radius for focus 2
pd      = p(3:end-vflag-expflag); % diffusion coefficients initial guesses
expfact = 1./p(numel(pd)+3:end-vflag); % triplet lifetime (s)
expfact = expfact(:).'; 
vel     = p(end-vflag+1:end); % velocity in um/s

t=t(:);
[m,~] = size(y);
ym =zeros(m,4);
ym(:,1)=y(:,1);
ym(:,2)=y(:,2);
ym(:,3)=(y(:,3)+y(:,5)); % summing up the correlations L1D1 X L2D2 and L1D2 X L2D1
ym(:,4)=(y(:,4)+y(:,6)); % summing up the correlations L2D1 X L1D2 and L2D2 X L1D1
% [m,n]=size(ym);
denom1 = 4*pd*t+w1^2; % denominator in the function of focus 1
denom2 = 4*pd*t+w2^2; % denominator in the function for focus 2
denomm = (denom1+denom2)/2; % denominator for cross correlation term
c=zeros(2*expflag+vflag,4);
% First laser auto
zz1 = [ones(m,1) (repmat(exp((-(vel.*t).^2)./denom1)./sqrt(denom1),1,2*(expflag))).*[ones(m,expflag) exp(-t*expfact)]];
c(:,1) = lsqnonneg(zz1,ym(:,1));
z(:,1) = zz1*c(:,1);
% Second laser auto
zz2 = [ones(m,1) (repmat(exp((-(vel.*t).^2)./denom2)./sqrt(denom2),1,2*(expflag))).*[ones(m,expflag) exp(-t*expfact)]];
c(:,2) = lsqnonneg(zz2,ym(:,2));
z(:,2) = zz2*c(:,2);
% First cross Second laser
zz3 = [ones(m,1) (repmat(exp((-(delta-vel.*t).^2)./denomm)./sqrt(denomm),1,2*(expflag))).*[ones(m,expflag) exp(-t*expfact)]];
c(:,3) = lsqnonneg(zz3,ym(:,3));
z(:,3) = zz3*c(:,3);
% Second cross First laser
zz4 = [ones(m,1) (repmat(exp((-(delta+vel.*t).^2)./denomm)./sqrt(denomm),1,2*(expflag))).*[ones(m,expflag) exp(-t*expfact)]];
c(:,4) = lsqnonneg(zz4,ym(:,4));
z(:,4) = zz4*c(:,4);

zz=[zz1;zz2;zz3;zz4];

if pic==1
    subplot(4,1,1:3)
    semilogx(t, ym, 'ob', t, z, 'r');
    axis tight
    subplot(4,1,4)
    semilogx(t, (ym-z)./sqrt(abs(z)),t,0*t)
    axis tight
    drawnow
elseif pic==2
    semilogx(t, ym, 'o', t, z); drawnow
else
    semilogx(t, ym, 'o', t, z); drawnow
end
% 
%  err = sum(sum((ym-z).^2./abs(z)));
%  err = sum((ym-z).^2);
ind = ym>0;
err = sum(ym(ind).*log(ym(ind)./z(ind))-ym(ind)+z(ind));



