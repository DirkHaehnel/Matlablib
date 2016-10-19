function matrixel = fnodipol(nm,lm,mm,n,l,m,pol,omega,x,h)
%
% Rectangle integration of <nm,lm,mm| exp(ik.r) \vec(pol).\grad |nlm >.
% Input variables:
% nm,lm,mm, n,l,m integers; Hydrogenic quantum numbers.
% pol;  is the polarisation vector which may have one of three values:
%        pol=0 for linear polarised light, (laser along the x-axis)
%        pol=-1 for right hand circular polarisation (x-iy) (laser along the z-axis) 
%        pol=+1 for left hand circular polarisation (x+iy)           
% 
% omega: frequency of light field. The wavenumber k is always taken to be k=omega/c
% x = vector of spatial mesh points, ie. x=0.001:h:1500;
% h: step size;
% WARNING: do never use x(1) =0  but always > 0, ie. 0.0001   
% Calling Routines:
% Routines in use: mlrange.m, spheint.m, radintR.m, radintddR.m

im = sqrt(-1);
c = 137.03587;
k= omega / c;
%k=0;
pol2=pol*pol;
if (pol2 ~= 0) &  (pol2 ~= 1)
 fprintf('I dont understand the value of polarisation, ='), pol
else 
  if pol == 0                 
% Linear polarized light (k=(kx,0,0))                
    vk1=pi/2;
    vk2=0;

    sum1 = 0; sum2 = 0;
    prod1=0; prod2=0;

    [lval,ml,nol]=  mlrange(lm,mm,l+1,m);
    for i=1:nol
      lam=lval(i); 
      mu=ml;
      spint1 = spheint(lm,mm,lam,mu,l+1,m);
      if spint1 ~=0
        prefac1 = 4*pi*im^lam*sqrt( ((l+m+1)*(l-m+1)) / ((2*l+3)*(2*l+1)) );
        ylmk = (-1)^mu*ylm(lam,mu,vk1,vk2);
        y=k*x; nu=lam+1/2;
        fac1=y.^(-0.5).*sqrt(pi/2).*besselj(nu,y);
        y1=x.*radial_funk(nm,lm,x).*radial_funk(n,l,x).*fac1; 
        y2=radial_funk(nm,lm,x).*radial_funkdiff(n,l,x).*fac1;
        intradR = sum(y1)*h; intradddR= sum(y2)*h;
        prod1 = prefac1*ylmk*spint1*(intradddR-l*intradR);
      end; 
      sum1=sum1+prod1;
    end;

    if l>0 
      [lval,ml,nol]=  mlrange(lm,mm,l-1,m);
      for i=1:nol
        lam=lval(i);
        mu=ml;
        spint2 = spheint(lm,mm,lam,mu,l-1,m);
        if spint2 ~=0
          prefac2 = 4*pi*im^lam*sqrt( ((l-m)*(l+m)) / ((2*l+1)*(2*l-1)) );
          ylmk = (-1)^mu*ylm(lam,mu,vk1,vk2);
          y=k*x; nu=lam+1/2;
          fac1=y.^(-0.5)*sqrt(pi/2).*besselj(nu,y);
          y1=x.*radial_funk(nm,lm,x).*radial_funk(n,l,x).*fac1; 
          y2=radial_funk(nm,lm,x).*radial_funkdiff(n,l,x).*fac1;
          intradR = sum(y1)*h; intradddR= sum(y2)*h;
          prod2 = prefac2*ylmk*spint2*(intradddR+(l+1)*intradR);
        end;  
        sum2=sum2+prod2;
      end;
    end;
  else             
% Circular polarized light (k=(0,0,kz)):              
    vk1=0;
    vk2=0;
    sum1 = 0; sum2 = 0;
    prod1=0; prod2=0;

    [lval,ml,nol]=  mlrange(lm,mm,l+1,m+pol);
    for i=1:nol
      lam=lval(i);
      mu=ml;
      spint1 = spheint(lm,mm,lam,mu,l+1,m+pol);
      if spint1 ~=0
        prefac1 = pol*4*pi*im^lam*sqrt( ((l+pol*m+1)*(l-pol*m+1)) / ((2*l+3)*(2*l+1)) );
        ylmk = (-1)^mu*ylm(lam,mu,vk1,vk2);
        y=k*x;
        nu=lam+1/2;
%        fac1=y.^(-0.5).*sqrt(pi/2).*besselj(nu,y);
        fac1=sphebess_func(lam,y);
        y1=x.*radial_funk(nm,lm,x).*radial_funk(n,l,x).*fac1;
        y2=radial_funk(nm,lm,x).*radial_funkdiff(n,l,x).*fac1;
        intradR = sum(y1)*h; intradddR= sum(y2)*h;
        prod1 = prefac1*ylmk*spint1*(intradddR-l*intradR);
      end;
      sum1=sum1+prod1;
    end;

    if l>0 
      [lval,ml,nol]=  mlrange(lm,mm,l-1,m+pol);
      for i=1:nol
        lam=lval(i);
        mu=ml;
        spint2 = spheint(lm,mm,lam,mu,l-1,m+pol);
        if spint2 ~=0
          prefac2 = -pol*4*pi*im^lam*sqrt( ((l-pol*m)*(l-pol*m-1)) / ((2*l+1)*(2*l-1)) );
          ylmk = (-1)^mu*ylm(lam,mu,vk1,vk2); 
          y=k*x; nu=lam+1/2;
%          fac1=y.^(-0.5)*sqrt(pi/2).*besselj(nu,y);
          fac1=sphebess_func(lam,y);
          y1=x.*radial_funk(nm,lm,x).*radial_funk(n,l,x).*fac1;
          y2=radial_funk(nm,lm,x).*radial_funkdiff(n,l,x).*fac1;
          intradR = sum(y1)*h; intradddR= sum(y2)*h;
          prod2 = prefac2*ylmk*spint2*(intradddR+(l+1)*intradR);
        end;
        sum2=sum2+prod2;
      end;
    end;
  end
  matrixel = sum1+sum2;
end
