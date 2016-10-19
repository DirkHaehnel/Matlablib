% read the IRF
fin = fopen('c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\PumpProbe\l07092490_IRF_P.txt','r');
for j=1:3 fgetl(fin); end
x = fscanf(fin,'%f %f\n',inf);
fclose(fin);
x = [x(1:2:end) x(2:2:end)];

% fin = fopen('c:\Joerg\Doc\Meixner\Frank\PumpProbe\l07092492_IRF_PPr.txt','r');
% for j=1:3 fgetl(fin); end
% y = fscanf(fin,'%f %f\n',inf);
% fclose(fin);
% y = [y(1:2:end) y(2:2:end)];

% only pump beam
cd c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\PumpProbe\2ns_RG500_P
fnames = dir('*.txt');
for j=1:length(fnames)
    fin = fopen(fnames(j).name,'r');
    for k=1:3 
        fgetl(fin); 
    end
    tmp = fscanf(fin,'%f %f\n',inf);
    fclose(fin);
    y(:,:,j) = [tmp(1:2:end) tmp(2:2:end)];
end    
    
% pump + probe beam
cd c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\PumpProbe\2ns_RG500_PPr
fnames = dir('*.txt');
for j=1:1
    fin = fopen(fnames(j).name,'r');
    for k=1:3 
        fgetl(fin); 
    end
    tmp = fscanf(fin,'%f %f\n',inf);
    fclose(fin);
    z(:,:,j) = [tmp(1:2:end) tmp(2:2:end)];
end

cd c:\Joerg\Doc\Nanophotonics&Optics\FrankSchleifenbaum\PumpProbe
ind = x(:,1)>35 & x(:,1)<60; % restrict analysis to times between 35 and 60 ns
dt = mean(diff(x(:,1))); % TCSPC channel width
p = 100/dt; % periodicity of TCSPC excitation
tp = (1:p)';
tau1 = [0.1 1]/dt; m = 2;
tau2 = [1]/dt; 
t = 1:sum(ind); n = length(t);
irf = x(ind,2);
fac = 0.5;
for j=1:1%size(y,3)
    close all
    param1 = [0.1/dt; tau1'];
    % Decay times and Offset are assumed to be positive.
    paramin = [-inf, zeros(1,length(tau1))];
    paramax = [inf, inf*ones(1,length(tau1))];
    param1 = Simplex('lsfit', param1, paramin, paramax, [], [], irf, y(ind,2,j), p);
    [err A1] = lsfit(param1, irf, y(ind,2,j), p);
    
    param2 = [3/dt; param1(2:end); tau2'];
    % Decay times and Offset are assumed to be positive.
    paramin = [-inf, param1(2:end)', zeros(1,length(tau2))];
    paramax = [inf, param1(2:end)', inf*ones(1,length(tau2))];
    param2 = Simplex('FrankPumpProbeFitFun', param2, paramin, paramax, [], [], irf, z(ind,2,j), p, param1(1), A1(2:end), param1(2:end));

    c1 = param1(1); res(j).c1 = dt*c1;
    c2 = param2(1); res(j).c2 = dt*c2;
    tau1 = param1(2:end)'; res(j).tau1 = dt*tau1;
    tau2 = param2(2:end)'; res(j).tau2 = dt*tau2;
    x1 = exp(-(tp-1)*(1./tau1))*diag(1./(1-exp(-p./tau1)));
    x2 = exp(-(tp-1)*(1./tau2))*diag(1./(1-exp(-p./tau2)));

    irs1 = (1-c1+floor(c1))*irf(rem(rem(t-floor(c1)-1, n)+n,n)+1) + (c1-floor(c1))*irf(rem(rem(t-ceil(c1)-1, n)+n,n)+1);
    irs2 = (1-c2+floor(c2))*irf(rem(rem(t-floor(c2)-1, n)+n,n)+1) + (c2-floor(c2))*irf(rem(rem(t-ceil(c2)-1, n)+n,n)+1);
    z1 = convol(irs1, x1);
    z2 = convol(irs2, x2);
    z2 = [ones(size(z1,1),1) z1 z2];
    z1 = [ones(size(z1,1),1) z1];
    A1 = z1\y(ind,2,j); res(j).A1 = A1;
    A2 = z2\z(ind,2,j); res(j).A2 = A2;
    %A = lsqnonneg(z,y);
    z1 = z1*A1;
    z2 = z2*A2;

    t = dt*t;
    tau1 = dt*tau1;
    tau2 = dt*tau2;
    subplot('position',[0.1 0.4 0.8 0.5])
    % 	semilogy(t, min(irs1) + (irs1-min(irs1))/(max(irs1)-min(irs1))*max(y(ind,2,j)), t, y(ind,2,j), 'o', t, z1, ...
    %         t, min(irs2) + (irs2-min(irs2))/(max(irs2)-min(irs2))*max(z(ind,2,j)), t, z(ind,2,j), 'o', t, z2);
 	semilogy(t, y(ind,2,j), 'o', t, z1, t, z(ind,2,j), 'o', t, z2);
        
    xlabel('Time in ns');
    ylabel('Log Count');
    v = axis;
    v(1) = min(t);
    v(2) = max(t);
    axis(v);

    s = sprintf('offset = %3.3f   %3.3f', z1(1,1), z2(1,1));
    text(max(t)/2,fac*v(4),s);
    s = ['AMP1 = '];
    A1 = A1(2:end)/sum(A1(2:end));
    for k=1:length(A1)
        s = [s sprintf('%1.3f',A1(k)) '   '];
    end
    text(max(t)/2,fac^2*v(4),s);
    s = ['AMP2 = '];
    A2 = A2((2+length(tau1)):end)/sum(A2((2+length(tau1)):end));
    for k=1:length(A2)
        s = [s sprintf('%1.3f',A2(k)) '   '];
    end
    text(max(t)/2,fac^3*v(4),s);
    s = ['TAU1 = '];
    for k=1:length(tau1)
        s = [s sprintf('%3.3f',tau1(k)) '   '];
    end
    text(max(t)/2,fac^4*v(4),s);    
    s = ['TAU2 = '];
    for k=1:length(tau2)
        s = [s sprintf('%3.3f',tau2(k)) '   '];
    end
    text(max(t)/2,fac^5*v(4),s);    
    
    subplot('position',[0.1 0.1 0.8 0.2])
    plot(t,(y(ind,2,j)-z1)./sqrt(abs(z1)),t,(z(ind,2,j)-z2)./sqrt(abs(z2)));
    chi1 = sum((y(ind,2,j)-z1).^2./abs(z1))/(n-m);
    chi2 = sum((z(ind,2,j)-z2).^2./abs(z2))/(n-m);
    v = axis;
    v(1) = min(t);
    v(2) = max(t);
    axis(v);
    xlabel('Time in ns');
    ylabel('Residues');
    s = sprintf('%3.3f   &   %3.3f', chi1, chi2);
    text(max(t)/2,v(4)-0.1*(v(4)-v(3)),['\chi^2 = ' s]);
    set(gcf,'units','normalized','position',[0.01 0.05 0.98 0.83])

    eval(['print -dpng -r300 ' fnames(j).name(1:end-4)])
end

save FitResults res

