name = 'C:\joerg\temp\video006.avi';
start = 200; 
para = aviinfo(name);

ar = 35:105;
br = 5:80;

a = zeros(para.NumFrames-start,1);
b = a;

warning off
for j=1:para.NumFrames-start;
    mov=aviread(name,start+j);
    x = double(mov.cdata(:,:,1));
    x = x(ar,br);
    x=max(max(x))-x;
    y = x>1.25*mconv2(x,disk(10));
    y = erode(y,11,5);
    a(j) = (1:size(y,1))*sum(y')'/sum(sum(y));
    b(j) = (1:size(y,2))*sum(y)'/sum(sum(y));
    y(a(j),b(j)) = inf;
    %mim(y); drawnow
    if mod(j,10)==0
        plot(a(1:j),b(1:j)); axis image; drawnow;
    end
end
warning on
    
break

n = 250;
clear ha; ht=-25:0.1:25; for j=1:n ha(:,j) = mhist(a(j+1:end)-a(1:end-j),ht); end
clear hb; ht=-25:0.1:25; for j=1:n hb(:,j) = mhist(b(j+1:end)-b(1:end-j),ht); end
ma = ht*ha./sum(ha); mb = ht*hb./sum(hb);
ga = sqrt(ht.^2*ha./sum(ha)-ma.^2);
gb = sqrt(ht.^2*hb./sum(hb)-mb.^2);

ka = simplex('OneMExpFun',1,0,[],[],[],1:n,ga.^2);
[err, ca, gat] = OneMExpFun(ka, 1:n, ga.^2);
kb = simplex('OneMExpFun',1,0,[],[],[],1:n,gb.^2);
[err, cb, gbt] = OneMExpFun(kb, 1:n, gb.^2);




