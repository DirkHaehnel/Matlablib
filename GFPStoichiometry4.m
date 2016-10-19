cd E:\Daten\Gregor\050428_GFP
names = dir('*ho400fcs.mat');

for j=1:1; %length(names)
    eval(['load ' names(j).name]);
    z = MFCS2Z(res.auto(:,:,:,1,2),res.intensity(:,1),res.intensity(:,2)) + ...
        MFCS2Z(res.auto(:,:,:,2,1),res.intensity(:,2),res.intensity(:,1));
    t = res.autotime;
    y(:,1) = z(:,1,1); 
    y(:,2) = z(:,1,2)+z(:,2,1);
    y(:,3) = z(:,2,2);
    y(:,4) = z(:,2,3)+z(:,3,2);
    
    close all
    for k=1:5e2
        tt = (1+k/1e2)*t;
        tt = tt(tt<max(t));
        tmp = interp1(t,y(:,1),tt,'cubic');
        c = tmp\y(1:length(tmp),2);
        semilogx(tt,y(1:length(tmp),2),tt,c*tmp,'o'); drawnow
        err(k) = sum((y(1:length(tmp),2)-c*tmp).^2);
    end
        

end
    

