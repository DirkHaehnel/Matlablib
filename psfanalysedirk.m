numofpix=40;
for i=1:numofpix
    name = ['\\jesrv\WG-Data\ISM\2012-03\26\MessungPSF\AreaISM_X-Y_' num2str(i) '_260312_1412_ISM.mat'];
    load(name);
    image=reshape(results.SHR,16,16,150,150);
    ergebniss1(i,:,:) = squeeze(sum(sum(image,1),2));
    ergebniss2(i,:,:) = squeeze(sum(sum(image,3),4));
    [mx my wx wy amp]=Gauss2D([],[],squeeze(ergebniss2(i,5:40,10:45)));
    abm=size(results.SHR);
    ergebniss(i,:) = [results.Info.zPos (mx+5)./abm(1) (my+10)./abm(2) wx*0.08 wy*0.08 amp];
end

subplot(2,3,1);
plot(ergebniss(:,1),[ergebniss(:,2) ergebniss(:,3)]);
legend({'mx','my'});
xlabel('z-Pos / [nm]');
ylabel('Position / [nm]');
subplot(2,3,2);
plot(ergebniss(:,1),[ergebniss(:,4) ergebniss(:,5)]);
legend({'wx','wy'});
xlabel('z-Pos / [nm]');
ylabel('1/e-Gauss, / [\mum]');
subplot(2,3,3);
plot(ergebniss(:,1),ergebniss(:,6));
legend({'amp'});
xlabel('z-Pos / [nm]');
ylabel('Counts');
drawnow
while 1==1
for i=1:numofpix
    subplot(2,3,4)
    mim(squeeze(ergebniss2(i,:,:)))
    subplot(2,3,5)
    mim(squeeze(ergebniss1(i,:,:)))
    drawnow
end
end