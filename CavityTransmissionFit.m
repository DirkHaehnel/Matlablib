function [err tw c y] = CavityTransmissionFit(d, lambda, trans, nwater, met, d0, d1)

if nargin<7
    d1 = d0;
    d0 = d(2);
    d = d(1);
end

nglass = 1.52;
if ischar(met)
    load metals
    %eval(['nag = interp1(wavelength,' met ',lambda);'])
else
    nag = interp1(met(:,1),met(:,2),lambda);
end

% transmission calculation
back = ones(numel(lambda),1);
tw = back;
for j=1:length(lambda)
    [~, ~, tw(j)] = Fresnel(nglass,[nglass nag(j) nwater nag(j) nglass],[d0 d d1]/lambda(j)*2*pi);
end
tw = abs(tw).^2;

c = [back lambda.^2/sum(lambda.^2) tw]\trans;

y = [back lambda.^2/sum(lambda.^2) tw]*c;
plot(lambda, trans, lambda, y); drawnow

err = sum((trans-y).^2);

return

dname = 'c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Converted transmissions\R6G in water 20111121\'
fnames = dir([dname '*.asc']);
load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\ModelSilver488ex30nm.mat','model'); 
ind = 5:25;
p = zeros(2,numel(fnames));
for j=1:length(fnames)
    x = load([dname fnames(j).name]);
    mm = x(x(:,2)==max(x(:,2)),1);
    close
    p(:,j) = Simplex('CavityTransmissionFit', [polyval(polyfit(model.tmax(ind),model.dv(ind),1),mean(mm)) 30], [], [], [], [], x(:,1), x(:,2), 1.33, 'silver', 85);
    num(j) = str2num(fnames(j).name(1:end-4));
    xlabel('wavelength (nm)'); ylabel('transmission (a.u.)'); 
    ax = axis;
    s = {['\itd\rm = ' mnum2str(p(1,j),3,1) ' nm'],['\itd\rm_{Ag} = ' mnum2str(p(2,j),3,1) ' nm']};
    text(ax(1)+0.7*diff(ax(1:2)),ax(3)+0.8*diff(ax(3:4)),s);
    eval(['print -dpng -r300 ' '''' dname fnames(j).name(1:end-4) ''''])
end

return 

dname = 'c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Converted transmissions\R6G in glycerol 20111121\'
fnames = dir([dname '*.asc']);
load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\ModelSilver488exGlycerol30nm.mat','model'); 
ind = 5:25;
p = zeros(2,numel(fnames));
for j=1:length(fnames)
    x = load([dname fnames(j).name]);
    mm = x(x(:,2)==max(x(:,2)),1);
    close
    p(:,j) = Simplex('CavityTransmissionFit', [polyval(polyfit(model.tmax(ind),model.dv(ind),1),mean(mm)) 30], [], [], [], [], x(:,1), x(:,2), 1.47, 'silver', 85);
    num(j) = str2num(fnames(j).name(1:end-4));
    xlabel('wavelength (nm)'); ylabel('transmission (a.u.)'); 
    ax = axis;
    s = {['\itd\rm = ' mnum2str(p(1,j),3,1) ' nm'],['\itd\rm_{Ag} = ' mnum2str(p(2,j),3,1) ' nm']};
    text(ax(1)+0.7*diff(ax(1:2)),ax(3)+0.8*diff(ax(3:4)),s);
    eval(['print -dpng -r300 ' '''' dname fnames(j).name(1:end-4) ''''])
end
