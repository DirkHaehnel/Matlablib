close all

clear c ctrs sumd T err
for j=1:size(amp1,2) 
    subplot(121);
    
    if 1
        c(:,j) = polyfit(amp1(:,j),amp2(:,j),1);
        err(j) = sqrt(sum((amp2(:,j)-polyval(c(:,j),amp1(:,j))).^2));
        [T,ctrs,sumd(j,:)] = kmeans([amp1(:,j), amp2(:,j)],2);
        dd(j) = sqrt(sum(diff(ctrs).^2));
        scatter(amp1(:,j),amp2(:,j),10,T,'filled');
        hold on
        plot(sort(amp1(:,j)),polyval(c(:,j),sort(amp1(:,j))));
        hold off
    end
    
    axis('equal');
    subplot(122);
    plot(1:size(amp1,1),[amp1(:,j) amp2(:,j)]); title(j);
    drawnow
    %pause; 
end

if 0
    [a,bin] = mHist(c(1,:),[-0.2:0.05:1.2]);
    para = [0 0.5 0.1 0.1];
    close;  for j=1:10 para=Simplex('Gauss',para,[],[],[],[],bin,a,0,[],1,exp(-a/10)); end
    [cc, cc] = Gauss(para,bin,a,0);
    t = -0.2:0.001:1.2;
    y = [exp(-(t-para(1)).^2/2/para(3)^2)' exp(-(t-para(2)).^2/2/para(4)^2)'];
    bar(bin,a);
    hold on
    plot(t,y*cc,'r');
    hold off
end

