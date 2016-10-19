function [auto, autotime] = normtttr2xfcs(y, num, Ncasc, Nsub)

% Function [auto, autotime] = tttr2xfcs(y, num, Ncasc, Nsub) calculates
% the correlation and crosscorrelation between photon streams with arrival 
% times y, where num is a matrix indicating which photon corresponds to
% which photon stream


y = round(y(:));
if size(num,1)<size(num,2)
    num = num';
end
autotime = cumsum(reshape(repmat(2.^(0:Ncasc-1),Nsub,1),Ncasc*Nsub,1));
dt = max(y)-min(y); % Total time of the experiment T.
auto = zeros(Ncasc*Nsub,size(num,2),size(num,2));
shift = 0;
delta = 1; % time resolution
ind=num~=0;
tphot=sum(ind,1); % Total number of photons in each detector
for j=1:Ncasc
    [y, k] = unique(y);
	tmp = cumsum(num);
    num = diff([zeros(1,size(num,2)); tmp(k,:)]);
    for k=1:Nsub
        shift = shift + delta;
        lag = round(shift/delta);
        [~,i1,i2] = cIntersect(y,lag);
        auto(k+(j-1)*Nsub,:,:) = num(i1,:)'*num(i2,:)/delta*dt/(dt-autotime(k+(j-1)*Nsub)); % <I1(t).I2(t+lag)>
    end
    y=round(0.5*y);
    delta = 2*delta;
end
% n_detectors = size(num,2);
% for i=1:n_detectors
%     for j=1:n_detectors
%         norm = (tphot(i).*tphot(j));
%         auto(:,i,j) = auto(:,i,j)*dt/norm;
%     end
% end
% 

norm=shiftdim(repmat((tphot.'*tphot),[1,1,numel(autotime)]),2); % <I1(t)>.<I2(t+lag)>*T
auto=auto*dt./norm; % <I1(t).I2(t+lag)>/<I1(t)>.<I2(t+lag)>