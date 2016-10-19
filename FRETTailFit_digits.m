function [fretE, res] = FRETTailFit(dt, donor_exp, fret_exp, tau, t1, t2)

set(0,'defaultFigurePosition',[100 100 660 520].*[0.2 0.5 1.5 1.2])
set(0,'defaultAxesFontName', 'Times', 'defaultAxesFontSize', 16, 'defaultLineLineWidth', 2, ...
    'defaultTextFontName', 'Times', 'defaultTextFontSize', 16);
set(0,'defaultaxeslinewidth',1.5)
set(0,'defaultAxesColorOrder', [1 0 0; 0 0 1; 0 0.5 0; 0.75 0.75 0; 0 0.75 0.75; 0.75 0 0.75; 0.25 0.25 0.25])
set(0,'defaultAxesLineWidth', 1.5)
set(0,'defaultLineMarkerSize',3)

if nargin<4 || isempty(tau)
    tau = [1 2];
end
if nargin<5 || isempty(t1)
    t1 = 0.5;
end
if nargin<6 || isempty(t2)
    t2 = 12;
end
mm = size(donor_exp,2);
if ~(size(fret_exp,2)==mm)
    print('Unequal number of donor and fret curves!')
    return
end
t = (round(t1/dt):round(t2/dt))*dt;
t = t-min(t);
[~, idx] = max(donor_exp);
for j=1:mm
    donor(:,j) = donor_exp((idx(j)+round(t1/dt)):(idx(j)+round(t2/dt)),j);
end
[~, idx] = max(fret_exp);
for j=1:mm
    fret(:,j) = fret_exp((idx(j)+round(t1/dt)):(idx(j)+round(t2/dt)),j);
end

irf = zeros(size(donor,1),1); 
irf(1)=1; 

for j=1:mm
    %[~, tau] = DistFluofit(irf, donor(:,j), max(t)/dt, dt, [0 0], 1, 1);
    pdonor(:,j) = Simplex('ExpFun',tau,[],[],[],[],t,donor(:,j),1);
    [~, cdonor(:,j), ~, zdonor(:,j)] = ExpFun(pdonor(:,j),t,donor(:,j),1);
    subplot(4,1,1:3)
    title(['donor curve ' int2str(j)])
    ylabel('fluorescence intensity (a.u.)')
    ax = axis;
    for k=1:size(pdonor,1)
        text(ax(1)+0.5*(ax(2)-ax(1)),ax(3)+(0.9-0.1*k)*(ax(4)-ax(3)),['\tau_ ' int2str(k) ' = ' mnum2str(pdonor(k,j),1,2) ' ns'])
    end
    subplot(4,1,4)
    xlabel('time (ns)')
    ylabel('residual')
    eval(['print -dpng -r300 DonorFit' int2str(j)])
    figure
    %[~, tau] = DistFluofit(irf, fret(:,j), max(t)/dt, dt, [0 0], 1, 1);
    pfret(:,j) = Simplex('ExpFun',tau,[],[],[],[],t,fret(:,j),1);
    [~, cfret(:,j), ~, zfret(:,j)] = ExpFun(pfret(:,j),t,fret(:,j),1);
    subplot(4,1,1:3)
    title(['FRET curve ' int2str(j)])
    ylabel('fluorescence intensity (a.u.)')
    ax = axis;
    for k=1:size(pfret,1)
        text(ax(1)+0.5*(ax(2)-ax(1)),ax(3)+(0.9-0.1*k)*(ax(4)-ax(3)),['\tau_ ' int2str(k) ' = ' mnum2str(pfret(k,j),1,2) ' ns'])
    end
    subplot(4,1,4)
    xlabel('time (ns)')
    ylabel('residual')
    eval(['print -dpng -r300 FRETFit' int2str(j)])
    figure 
end
fretE = (sum(cfret(2:end,:))./sum(cfret(2:end,:).*pfret(:,:))-sum(cdonor(2:end,:))./sum(cdonor(2:end,:).*pdonor(:,:)))./(sum(cfret(2:end,:))./sum(cfret(2:end,:).*pfret(:,:)));
plot(t,(mean(fret,2)-min(mean(zfret,2)))/(max(mean(zfret,2))-min(mean(zfret,2))),'o',...
    t,(mean(donor,2)-min(mean(zdonor,2)))/(max(mean(zdonor,2))-min(mean(zdonor,2))),'o',...
    t,(mean(zfret,2)-min(mean(zfret,2)))/(max(mean(zfret,2))-min(mean(zfret,2))),'r',...
    t,(mean(zdonor,2)-min(mean(zdonor,2)))/(max(mean(zdonor,2))-min(mean(zdonor,2))),'b','MarkerSize',3);
legend({'FRET','donor only'})
xlabel('time (ns)')
ylabel('fluorescence intensity (normalized)')
text(0.5*(t2-t1),0.7,['FRET efficiency = ' mnum2str(mean(fretE),1,3)])
axis([0 t2-t1 0 1.1])
print -dpng -r300 FRETEfficiency

res.t = t; % time axis
res.donor = donor; % donor curve tail
res.zdonor = zdonor; % fitted donor curve
res.pdonor = pdonor; % fitted donor lifetime(s)
res.cdonor = cdonor; % amplitude(s) of donor lifetime(s)
res.fret = fret; % fret curve tail
res.zfret = zfret; % fitted fret curve
res.pfret = pfret; % fret lifetime(s)
res.cfret = cfret; % amplitude(s) of fret lifetime(s) 
res.fretE = fretE; % FRET efficienies

% Write curves into ASCII file

tmp = [t(:) donor zdonor fret zfret];
save('data.txt', 'tmp', '-ASCII');

