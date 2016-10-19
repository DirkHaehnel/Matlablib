% FCSFit

close all
AvogadroConstant = 6.0221419947e23;

path = 'D:\Joerg\Doc\FCS\2Focus\CAM\2006-04-24\'
names = dir([path 'CAM_WT*.mat']);
ind = [];
for j=1:length(names)
    if ~isempty(findstr(names(j).name,'_parts'))
        ind = [ind j];
    end
end
%names(ind) = [];
names = names(ind);

% names = dir([path 'CAMT27C*.mat'])

for j=1:length(names)
    pos = findstr(names(j).name,'P');
    zd(j) = str2num(names(j).name(pos+1:pos+2));
end
[zd ind] = sort(zd);
names = names(ind);

% CaMconc = [12, 62, 113, 205, 285, 393, 566, 830, 1420, 2770, 240000, 1e6 10e3 100e3 200e3];
% zd = CaMconc(zd);

av = 100/58;
lam = [0.64 0.67]/1.33;
dist = 0.4;

global pd
expflag = 1;
clear p v dc trip conc
for j=1:length(names)
    [auto,t] = FCSCrossRead([path names(j).name],[1e-6 5e-1]);
%     eval(['load ' path names(j).name ' parameters']); 
%     if isstr(parameters.Temperature)
%         TT(j) = str2num(parameters.Temperature);
%     else
%         TT(j) = parameters.Temperature;
%     end

    if size(auto,3)>1
        for k=1:size(auto,3)
            pd = [1 1e-3 1e-3 1e-3];
            p(:,j,k) = simplex('GaussFCS',[0.37 0.12],zeros(1,2),[],[],[],av,lam,dist,1e2*t,auto(:,1:2,k),auto(:,3,k),expflag);
            [err, pd, c] = GaussFcs(p(:,j,k),av,lam,dist,1e2*t,auto(:,1:2,k),auto(:,3,k),expflag);

            v(j,k) = GaussDetectionVolume(p(:,j,k),av,lam);
            dc(j,k) = 1./pd(1)*1e-6
            if ~isempty(expflag)
                trip(:,j,k) = pd(2:end)*1e-2;
            end
            conc(j,k,:) = c(2,:)./c(1,:)/1e-15/AvogadroConstant;
        end
    else
        pd = [1 1e-3 1e-3 1e-3];
        p(:,j) = simplex('GaussFCS',[0.37 0.12],zeros(1,2),[],[],[],av,lam,dist,1e2*t,auto(:,1:2),auto(:,3),expflag);
        [err, pd, c] = GaussFcs(p(:,j),av,lam,dist,1e2*t,auto(:,1:2),auto(:,3),expflag);

        v(j) = GaussDetectionVolume(p(:,j),av,lam);
        dc(j) = 1./pd(1)*1e-6
        if ~isempty(expflag)
            trip(:,j) = pd(2:end)*1e-2;
        end
        conc(j,:) = c(2,:)./c(1,:)/1e-15/AvogadroConstant;

    end

end

load viscosity
vis = interp1(T,visc(1,:),273.15+TT);
vis0 = interp1(T,visc(1,:),273.15+25);
dcc = dc.*vis./(TT+273.15)*(25+273.15)/vis0

save CaMFCSFit zd dc p names av lam dist dc dcc trip v conc TT err

return

for j=1:length(names)
    [auto,t] = FCSCrossRead([path names(j).name],[1e-6 5e-1]);

    pd = [1 1e-3 1e-3 1e-3];
    p(:,j) = simplex('GaussFCS',p(:,j),zeros(1,2),[],[],[],av,lam,dist,1e2*t,auto(:,1:2),auto(:,3),expflag);

    [err, pd, c] = GaussFcs(p(:,j),av,lam,dist,1e2*t,auto(:,1:2),auto(:,3),expflag);

    v(j) = GaussDetectionVolume(p(:,j),av,lam);
    dc(j) = 1./pd(1)*1e-6
    if ~isempty(expflag)
        trip(:,j) = pd(2:end)*1e-2;
    end
    conc(j,:) = c(2,:)./c(1,:)/1e-15/AvogadroConstant;
    
end


return

conc0=10.^(-9:0.1:-2); 
ind=[1 0 0 1 0 1 0 1 0 1 0 1 1 0];
ind1=ones(1,14)==1;
ind1([2 14])=false;
[err,c,z1]=MichaelisMenten(p,zd(ind)*1e-9,dc(ind),conc0);
[err,c,z2]=MichaelisMenten(p,zd(~ind & ind1)*1e-9,dc(~ind & ind1),conc0);
semilogx(zd(ind)*1e-9,dc(ind),'o',zd(~ind & ind1)*1e-9,dc(~ind & ind1),'v',zd(~ind & ~ind1)*1e-9,dc(~ind & ~ind1),'d',conc0,z1,'r',conc0,z2,'b')
xlabel('Ca^{2+} concentration [M]'); ylabel('diffusion coefficient [cm^2/s]')
