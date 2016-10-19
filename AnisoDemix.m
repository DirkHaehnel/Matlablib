% demix = 0:0.01:pi/4;
% k = 1;
% clear tmp1 tmp2 tmp3
% for j=1:length(demix)
%     c = cos(demix(j));
%     s = sin(demix(j));
%     tmp1(j) = (1+2*c*s/3)^2;
%     tmp2(j) = 2*(1+2*c*s/3)*(1/3+2*c*s);
%     tmp3(j) = (1/3+2*c*s)^2;
% end
% 
% amp = [0.9565 0.2442];
% ind = 1:length(demix);
% ind = round(mean([ind(sum(tmp2./tmp1<=amp(1))),ind(sum(tmp3./tmp1<=amp(2)))]));
% plot(demix/pi*180,tmp2./tmp1,demix/pi*180,tmp3./tmp1,...
%     demix/pi*180,amp(1)*ones(size(demix)),':',...
%     demix/pi*180,amp(2)*ones(size(demix)),':',...
%     demix(ind)/pi*180*ones(1,101),(0:0.01:1)*max(tmp2./tmp1),':')
% axis tight
% xlabel('mixing angle (°)');
% ylabel('rel. amplitudes of correlation peaks')

load ConfoAnisoSat1

demix = 0:0.005:pi/4;
k = 1;
clear tmp1 tmp2 tmp3
for j=1:length(demix)
    c = cos(demix(j));
    s = sin(demix(j));
    tmp1(j) = (c^4*aniso(k).ux(1,1)+c^2*s^2*(aniso(k).uxr(1,1)+aniso(k).uy(1,1))+s^4*aniso(k).uyr(1,1))*...
        (c^4*aniso(k).uyr(1,1)+c^2*s^2*(aniso(k).uxr(1,1)+aniso(k).uy(1,1))+s^4*aniso(k).ux(1,1));
    tmp2(j) = (c^4*aniso(k).ux(1,1)+c^2*s^2*(aniso(k).uxr(1,1)+aniso(k).uy(1,1))+s^4*aniso(k).uyr(1,1))*...
        (c^4*aniso(k).uy(1,1)+c^2*s^2*(aniso(k).ux(1,1)+aniso(k).uyr(1,1))+s^4*aniso(k).uxr(1,1)) + ...
        (c^4*aniso(k).uxr(1,1)+c^2*s^2*(aniso(k).uyr(1,1)+aniso(k).ux(1,1))+s^4*aniso(k).uy(1,1))*...
        (c^4*aniso(k).uyr(1,1)+c^2*s^2*(aniso(k).uxr(1,1)+aniso(k).uy(1,1))+s^4*aniso(k).ux(1,1));
    tmp3(j) = (c^4*aniso(k).uy(1,1)+c^2*s^2*(aniso(k).uyr(1,1)+aniso(k).ux(1,1))+s^4*aniso(k).uxr(1,1))*...
        (c^4*aniso(k).uxr(1,1)+c^2*s^2*(aniso(k).ux(1,1)+aniso(k).uyr(1,1))+s^4*aniso(k).uy(1,1));
end

amp = [0.9565 0.2442];
ind = 1:length(demix);
ind = round(mean([ind(sum(tmp2./tmp1<=amp(1))),ind(sum(tmp3./tmp1<=amp(2)))]));
plot(demix/pi*180,tmp2./tmp1,demix/pi*180,tmp3./tmp1,...
    demix/pi*180,amp(1)*ones(size(demix)),':',...
    demix/pi*180,amp(2)*ones(size(demix)),':',...
    demix(ind)/pi*180*ones(1,101),(0:0.01:1)*max(tmp2./tmp1),':')
axis tight
xlabel('mixing angle (°)');
ylabel('rel. amplitudes of correlation peaks')

return

%lifetime effect

a0 = (aniso(k).ux(1,1))*(aniso(k).uyr(1,1));
b0 = (aniso(k).ux(1,1))*(aniso(k).uy(1,1)) + ...
        (aniso(k).uxr(1,1))*(aniso(k).uyr(1,1));
c0 = (aniso(k).uy(1,1))*(aniso(k).uxr(1,1));

a = (aniso(k).ux(1,1)+exp(-2)*aniso(k).uxr(1,1))*(aniso(k).uyr(1,1)+exp(-1)*aniso(k).uy(1,1));
b = (aniso(k).ux(1,1)+exp(-2)*aniso(k).uxr(1,1))*(aniso(k).uy(1,1)+exp(-2)*aniso(k).uyr(1,1)) + ...
        (aniso(k).uxr(1,1)+exp(-2)*aniso(k).ux(1,1))*(aniso(k).uyr(1,1)+exp(-1)*aniso(k).uy(1,1));
c = (aniso(k).uy(1,1)+exp(-2)*aniso(k).uyr(1,1))*(aniso(k).uxr(1,1)+exp(-1)*aniso(k).ux(1,1));

% calculate ACF from TCSPC

for j=1:length(tau) tst(j)=sum(tcspc(:,1).*tcspc(mod(j:1:end+j-1,end)+1,2)); end
bck = mean(tcspc(end-100:end,:));
for j=1:length(tau) tst1(j)=sum((tcspc(:,1)-bck(1)).*(tcspc(mod(j:1:end+j-1,end)+1,2)-bck(2))); end

