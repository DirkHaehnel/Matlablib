nsolvent = 1.33; % ref. index solvent
dagtop = 85; % top silver thickness
dagbottom = 35; % bottom silver thickness
lamex = 488; % excitaion wavelength

% ideal model calculation
dv = 100:5:250;
model  = AlexeyModel('silver', nsolvent, lamex, dagbottom, dagtop, dv); 
eval(['save ModelSilver488ex30nm.mat' mint2str(j,2) 'nm model']);

% fitting the transmission curves
dname = 'c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Converted transmissions\R6G in water 20111121\'; % folder with transmission curves
fnames = dir([dname '*.asc']); 
irf = load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\TransmissionGauge.txt'); % spectrometer calibration curve
load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\ModelSilver488ex30nm.mat','model'); % ideal model
load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\ag_ellipsometry.mat'); % ag_lam = wavelength, ag_re = real part, ag_im = imaginary part
ind = 5:25;
p = zeros(2,numel(fnames));
num = zeros(numel(fnames),1);
for j=1:length(fnames)
    x = load([dname fnames(j).name]);
    x(:,2) = (x(:,2)-min(x(:,2)))./interp1(irf(:,1),irf(:,2),x(:,1));
    mm = x(x(:,2)==max(x(:,2)),1);
    close
    p(:,j) = Simplex('CavityTransmissionFit', [polyval(polyfit(model.tmax(ind),model.dv(ind),1),mean(mm)) 30], [], [], [], [], x(:,1), x(:,2), nsolvent, [ag_lam ag_re+1i*ag_im], dagtop);
    num(j) = str2double(fnames(j).name(1:end-4));
    xlabel('wavelength (nm)'); ylabel('transmission (a.u.)');
    ax = axis;
    s = {['\itd\rm = ' mnum2str(p(1,j),3,1) ' nm'],['\itd\rm_{Ag} = ' mnum2str(p(2,j),3,1) ' nm']};
    text(ax(1)+0.7*diff(ax(1:2)),ax(3)+0.8*diff(ax(3:4)),s);
    eval(['print -dpng -r300 ' '''' dname fnames(j).name(1:end-4) '''' 'New']) % saving figure on hard disk
end

x = load([dname 'transmissions and decays.txt']); % relation between measurement numbers and max transmission wavelength; #1, #2, max tranmsission wavelength
dv = zeros(size(x,1),2); dag = dv; tmax = zeros(size(x,1),1);
for j=1:size(x,1)
    dv(j,:) = [p(1,num==x(j,1)) p(1,num==x(j,2))];
    dag(j,:) = [p(2,num==x(j,1)) p(2,num==x(j,2))];
    tmax(j)= x(j,3);
end
eval(['save ' '''' dname 'TransmissionFitNew' '''' ' tmax dv dag'])

% modeling
eval(['load ' '''' dname 'TransmissionFitNew' ''''])
dv = mean(dv,2);
dag = mean(dag,2);
model  = AlexeyModel([ag_lam ag_re+1i*ag_im], nsolvent, lamex, dag, dagtop, dv);
model.tmax = tmax;
eval(['save ' '''' dname 'modelNew' '''' ' model'])

% fitting the lifetime data
ltdata = load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Water R6G.txt'); % lifetime data: max. transmission wavelength, lifetime
spectrum = load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\spectrum R6G in water.txt'); % free dye emission spectrum
load('c:\Joerg\Doc\Nanophotonics&Optics\AlexeyChizhik\RotoDiffCavity\Converted transmissions\R6G in water 20111121\modelNew.mat','model');
ind = ltdata(:,1)>=510 & ltdata(:,1)<=650; % cutting off edges 
[paraw taufitw lamw tauw tau0w] = AlexeyCavityFitFull(ltdata(ind,:), spectrum, model, 1000); % paraw(1) = quantum yield, 1/6/paraw(2) = rotational diffusion time
lam = linspace(lamw(1),lamw(end));
plot(lamw,taufitw,lamw,tauw,'ob',lam,tau0w*ones(size(lam)),':b')
axis([505 655 1 4.5])
xlabel('max. transmission wavelength (nm)')
ylabel('lifetime (ns)')




