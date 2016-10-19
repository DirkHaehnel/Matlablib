subseqlength=15; %[s]

horder=4;
[fn,pn]=uigetfile('*.tif','MultiSelect','on');
if ~iscell(fn)
    tmp=cell(1,1);
    tmp{1}=fn;
    fn=tmp;
end
pbo=statusbar('Analyze sequences...');
n=1;
while n<=length(fn)
    stack=[];
    fileinfo = imfinfo([pn fn{n}]);
    Nframes=length(fileinfo);
    stack=zeros(fileinfo(1).Height,fileinfo(1).Width,Nframes);
    
    fig=statusbar('Read stack...');
    for m=1:Nframes
        if n*m==1
            stack=double(imread([pn fn{n}],m,'info',fileinfo));
        else
            stack(:,:,m)=double(imread([pn fn{n}],m,'info',fileinfo));
        end
        fig=statusbar(m/Nframes,fig);
        if isempty(fig)
            break;
        end
    end
    delete(fig);
    
    if (n+1<=length(fn)) & (length(strfind(fn{n+1},'X2.tif'))>0)
        fileinfo = imfinfo([pn fn{n+1}]);
        Nframes=length(fileinfo);
        
        fig=statusbar('Read stack, part 2...');
        for m=1:Nframes
            stack(:,:,end+1)=double(imread([pn fn{n+1}],m,'info',fileinfo));
            fig=statusbar(m/Nframes,fig);
            if isempty(fig)
                break;
            end
        end
        delete(fig);
        n=n+1;
    end
    

    mtrace=squeeze(mean(mean(stack)));
    figure; plot(mtrace);
    
    xlabel('Frames');
    ylabel('Brightness, AU')
    n=n+1;
    
%      close all;
 mtrace=squeeze(mean(mean(stack)));
 x=(1:length(mtrace));
 x=x';
 y=mtrace; 
%masking, background subtraction
 Stack=mean(stack,3);
 thrhold=median(Stack(:))+std(Stack(:));
 mask = Stack > thrhold;
%  imagesp(mask);
 mim(mask);
 Stack=stack; % modified mask display
 Stack=permute(Stack,[3 1 2]);
 Stack=Stack(:,mask);
 size(stack)
 sum(mask(:))
 a=size(Stack);
 Stack=mean(Stack,2);
 figure;
 line(x,Stack,'Color','green');
 background=(y(:)*31*101 - Stack*sum(mask(:)))/(31*101-sum(mask(:)));
 plot(x,Stack,'r');
 line(x,background,'Color','green');
 title('signal&background');
 figure;plot(x,Stack-background,'r');
 title('signal with subtracted background');
 y=Stack - conv(background,ones(10,1)/10,'same');
 figure;plot(x,y);
 title('signal with subtracted background convolved');
 Y=conv(y,ones(10,1)/10,'same');
 figure;plot(x,Y,'r'); 
 title('previous; convolved');
 
%exponential decay approximation and subtraction
 
% figure; plot(Y);
% model=inline('p(1)*exp(-x/p(2))+p(3)-Y','p','x','Y');
% % x=x(50:end-2000);
% % Y=Y(50:end-2000);
% param=lsqnonlin(model,[mean(mtrace) length(mtrace)/3 mean(mtrace)],[],[],[],x,Y);
% Y1 = model(param,x,0);
% line(x,Y1,'Color','red');
% resid=Y-Y1;
% figure;
% plot(x,resid,'r');
% tr=mean(resid); % crucial
% 
% plot(x,resid,'r',x,resid > tr,'k');
% 
% low=resid < tr;
% Plo=lsqnonlin(model,[mean(mtrace) length(mtrace)/3 mean(mtrace)],[],[],[],x(low),Y(low));
% figure; 
% plot(x,Y,'r',x,model(Plo,x,0),'k');
% Phi=lsqnonlin(model,[mean(mtrace) length(mtrace)/3 mean(mtrace)],[],[],[],x(~low),Y(~low));
% line(x,model(Phi,x,0));
% 
% figure;
% plot(x,Y-(model(Phi,x,0)+model(Plo,x,0))/2,'r');
% line(x,Y-Y1);
% Y=Y-(model(Phi,x,0)+model(Plo,x,0))/2;

%end of bleaching

%on-off points
 mask1=Y > mean(Y); %maybe manual insertion needed
%   mask1=Y > 170;
 
 line(x,mask1);
 title('on-off');
 
 mask1=diff(mask1);
 figure; 
 plot(mask1);
 on=find(mask1 > 0)
 off=find(mask1 < 0)
 
 a=size(on);
 if a(1)~=1
 j=1
    while j< a(1)-1
     if on(j+1) - on(j)<100
         z=find(off >= on(j) & off <= on(j+1));
         off=cat(1,off(1:z-1),off(z+1:end));
         
         on(j)= (on(j+1) + on(j))/2;
         a(1)=a(1)-1;
         on= cat(1, on(1:j), on(j+2:end));
        
     else    
     j=j+1  
     end
    end
 end
 
 a=size(off);
 if a(1)~=1
 j=1
 while j< a(1)-1
     if off(j+1) - off(j)<100
         z=find(on >= off(j) & on <= off(j+1));
         on=cat(1,on(1:z-1),on(z+1:end));
         
         off(j)= (off(j+1) + off(j))/2;
         a(1)=a(1)-1;
         
         off= cat(1, off(1:j), off(j+2:end));
         
     else
     j=j+1;
     end
 end
 end
 
% if smth wrong in the beginning of sequence
        on=on(2:end);
        off=off(2:end);

 nbursts=3;% manual insertion required 
 
 %approximation
 modelOn=inline('p(1)+p(2)*(1-exp(-max(0,x-p(3))/p(4)))-y','p','x','y');
 modelOff=inline('p(1)+p(2)*(1-exp(-max(0,x-p(3))/p(4)))-y','p','x','y');
%modelOff=inline('p(1)+p(2)*(1-exp(-(x-p(3))/p(4)))-y','p','x','y');
 figure; plot(x,Y,'r');
 Pon=[];
 Poff=[];
 

 area1=1200;
 area2=1500;
%  area3=500;
%  area4=1000;
 
 if on(1)<off(1)
    for j=1:nbursts
      Pon=cat(1, Pon, lsqnonlin(modelOn,[mean(Y(on:off)) mean(Y(off:on(2:end))) on(j) 5],[],[],[],x(on(j)-area1:on(j)+area1),Y(on(j)-area1:on(j)+area1)));
%     Pon=cat(1, Pon, lsqnonlin(modelOn,[-355 355 on(j) 1],[],[],[],x(on(j)-area:on(j)+area),Y(on(j)-area:on(j)+area)));
      Poff=cat(1, Poff, lsqnonlin(modelOff,[mean(Y(on:off)) mean(Y(off:on(2:end))) off(j) 5],[],[],[],x(off(j)-area2:off(j)+area2),Y(off(j)-area2:off(j)+area2)));
%     Poff=cat(1, Poff, lsqnonlin(modelOff,[1400 650 off(j) 5],[],[],[],x(off(j)-area:off(j)+area),Y(off(j)-area:off(j)+area)));
    end
 else
     for j=1:nbursts
       Pon=cat(1, Pon, lsqnonlin(modelOn,[mean(Y(on:off(2:end))) mean(Y(off:on)) on(j) 10],[],[],[],x(on(j)-area1:on(j)+area1),Y(on(j)-area1:on(j)+area1)));
       Poff=cat(1, Poff, lsqnonlin(modelOff,[mean(Y(on:off(2:end))) mean(Y(off:on)) off(j) 10],[],[],[],x(off(j)-area2:off(j)+area2),Y(off(j)-area2:off(j)+area2)));
%      Poff=cat(1, Poff, lsqnonlin(modelOff,[1400 650 off(j) 5],[],[],[],x(off(j)-area:off(j)+area),Y(off(j)-area:off(j)+area)));
     end
     
 end
 
 Pon
 Poff
 for j=1:nbursts
 line(x(on(j)-area1:on(j)+area1),modelOn(Pon(j,:),x(on(j)-area1:on(j)+area1),0));
 line(x(off(j)-area2:off(j)+area2),modelOff(Poff(j,:),x(off(j)-area2:off(j)+area2),0),'color', 'black');
 end
 
%  title('Fit. Curves');
    xlabel('Frames');
    ylabel('Brightness, AU')
 timeOn=mean(Pon(:,4));
 timeOff=mean(Poff(:,4));
 Onstd=std(Pon(:,4));
 Offstd=std(Poff(:,4));
 lowerlevel=mean(Pon(:,1));
 upperlevel=mean(Poff(:,1));
 
 if nbursts == 1
     Onstd=NaN;
     Offstd=NaN;
 end

 save([pn 'Analysis\'  '/' num2str(fn{n}(1:end-4)) '.mat'],'timeOn','timeOff','Onstd', 'Offstd', 'lowerlevel','upperlevel','mask','mtrace','x','Y', 'on', 'off', 'modelOn', 'modelOff','nbursts', 'Pon', 'Poff');
    
    
    
    pbo=statusbar(n/numel(fn),pbo);
    if isempty(pbo)
        break;
    end
end
delete(pbo);