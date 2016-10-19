% BeadScanRead

clear mx my wx wy amp dd

path = 'D:\Doc\Anastasia\Scanning\'
%path = 'C:\Temp\Scanning\scanning020507\';

names = dir([path 'atto655down2_*.t3r']);
%names = dir([path 'beads_*.t3r']);

read_fun = 'PI_Read5';
%read_fun = 'SCX200Read';
%read_fun = 'PIE710Read';

close all hidden
pix = 50;
for j=1:length(names)
    zd(j) = str2num(names(j).name(14:18));
    %zd(j) = str2num(names(j).name(7:10));
    %zd = 0;
end
[zd ind] = sort(zd);
names = names(ind);
for j=1:length(names)
    close all
    fname = [path names(j).name];

    eval(['[tag, tim, tch, bin] = ' read_fun '(fname,''2focus'');'])

%     if j==1
%         semilogy(bin,tch);
%         pos = ginput(2);
%     end

%     filter = bin>pos(1,1) & bin<pos(2,1);
%     eval(['[tag1, tim1] = ' read_fun '(fname,[filter'' filter'']);'])
%     filter = bin<pos(1,1) | bin>pos(2,1);
%     eval(['[tag2, tim2] = ' read_fun '(fname,[filter'' filter'']);'])

if j==1
%         im1(:,:,j) = sum(tim1,3);
%         im2(:,:,j) = sum(tim2,3);
    im1(:,:,j) = tag(:,:,1);
    im2(:,:,j) = tag(:,:,2);
    else
%         mx = size(tim1,1);
%         my = size(tim1,2);
%         im1(1:mx,1:my,j) = sum(tim1,3);
%         im2(1:mx,1:my,j) = sum(tim2,3);
        mx = size(tag,1);
        my = size(tag,2);
        im1(1:mx,1:my,j) = tag(:,:,1);
        im2(1:mx,1:my,j) = tag(:,:,2);
    end
    
    mm = 20;
    for k1=-mm:mm
        for k2=-mm:mm
            tst(mm+1+k1,mm+1+k2) = sum(sum(im1(mm+1+k1:end-mm+k1,mm+1+k2:end-mm+k2,j).*im2(mm+1:end-mm,mm+1:end-mm,j)));
        end;
    end

    [x,y] = meshgrid(-mm:mm,-mm:mm);
    xm = x(tst==max(tst(:)));
    ym = y(tst==max(tst(:)));
    [px(j), py(j), wx(j), wy(j), amp(j)] = Gauss2D(x(1,mm+1+(max(-mm,xm-10):min(mm,xm+10))),y(mm+1+(max(-mm,ym-10):min(mm,ym+10)),1),...
        tst(mm+1+(max(-mm,ym-10):min(mm,ym+10)),mm+1+(max(-mm,xm-10):min(mm,xm+10))));

end

dd = sqrt(px.^2+py.^2)*pix;


return

clear rr px py tst
mm = 20;
for s1=1:2
    for s2=1:2
        for j=1:size(im1,3)
            imm1 = im1((s1-1)*end/2+1:s1*end/2,(s2-1)*end/2+1:s2*end/2,:);
            imm2 = im2((s1-1)*end/2+1:s1*end/2,(s2-1)*end/2+1:s2*end/2,:);
            for k1=-mm:mm
                for k2=-mm:mm
                    tst(mm+1+k1,mm+1+k2) = sum(sum(imm1(mm+1+k1:end-mm+k1,mm+1+k2:end-mm+k2,j).*imm2(mm+1:end-mm,mm+1:end-mm,j)));
                end;
            end

            [x,y] = meshgrid(-mm:mm,-mm:mm);
            xm = x(tst==max(tst(:)));
            ym = y(tst==max(tst(:)));
            [px(j), py(j)] = Gauss2D(x(1,mm+1+(max(-mm,xm-10):min(mm,xm+10))),y(mm+1+(max(-mm,ym-10):min(mm,ym+10)),1),...
                tst(mm+1+(max(-mm,ym-10):min(mm,ym+10)),mm+1+(max(-mm,xm-10):min(mm,xm+10))));
        end
        rr(:,s1,s2) = sqrt(px.^2+py.^2)'*pix;
    end
end


return

mim(im1(:,:,1))
[b,a]=ginput(2);
a = round(a);
b = round(b);
for j=1:size(im1,3)
    [px1(j), py1(j), wx1(j), wy1(j), amp1(j)] = Gauss2D(im1(a(1):a(2),b(1):b(2),j));
    [px2(j), py2(j), wx2(j), wy2(j), amp2(j)] = Gauss2D(im2(a(1):a(2),b(1):b(2),j));    
end
rr=sqrt((px1-px2).^2+(py1-py2).^2)*pix
