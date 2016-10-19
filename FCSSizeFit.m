function [rad w0 a0 triplet c err] = FCSSizeFit(name,p,Tmeas,expflag,groupflag,para,pmin,pmax)

close all
global pd

if nargin<4 || isempty(expflag)
    expflag = 0;
end
if nargin>5 && ~isempty(para)
    av = para(1);
    lam = para(2:3);
    dist = para(4);
else
    av = 100e3/58;
    lam = [640 670]/1.33;
    dist = 403;
end

if isstruct(name)
    y = name.y;
    t = name.t;
else
    [y,t] = FCSCrossRead(name);
end
y(:,3,:) = (y(:,3,:)+y(:,4,:))/2;
y(:,4,:) = [];

if size(y,3)>1
    y = y(:,:,~isnan(sum(sum(y,1),2)));
end

if nargin<2 || isempty(p)
    p = [1 400 150]; % radius, (shell), beam waist, aperture parameter; 
end
if nargin<8 || isempty(pmin)
    pmin = 0*p;
end
if nargin<9 || isempty(pmax)
    pmax = inf*ones(size(p));
end
if nargin>4 && ~isempty(groupflag) && size(y,3)>1
    %tst = mean(abs(squeeze(sum(y))-mean(sum(y),3)'*ones(1,size(y,3)))./(mean(sum(y),3)'*ones(1,size(y,3))))<0.2;
    tst = ones(1,size(y,3))==1;
    p0 = Simplex('SizeFCS',p,pmin,pmax,[],[],av,lam,dist,Tmeas,1e6*t,sum(y(:,1:2,tst),3),sum(y(:,3:end,tst),3),expflag);
    if exist('pd')==1 pd0 = pd; end
    for j=1:ceil(size(y,3)/groupflag)
        if exist('pd')==1 pd = pd0; end        
        p = Simplex('SizeFCS',p0,pmin,pmax,[],[],av,lam,dist,Tmeas,1e6*t,sum(y(:,1:2,(j-1)*groupflag+1:min([j*groupflag,end])),3),sum(y(:,3:end,(j-1)*groupflag+1:min([j*groupflag,end])),3),expflag);
        [err(j), c(:,:,j)] = SizeFCS(p,av,lam,dist,Tmeas,1e6*t,sum(y(:,1:2,(j-1)*groupflag+1:min([j*groupflag,end])),3),sum(y(:,3:end,(j-1)*groupflag+1:min([j*groupflag,end])),3),expflag);
        if expflag>0
            triplet(:,j) = p(end-expflag+1:end);
        else 
            triplet = [];
        end
        rad(:,j) = p(1:2);
        w0(j) = p(3);
        a0(j) = p(4);
    end
    subplot(4,1,1:3)
    [tmp, tmp, z] = SizeFCS(p0,av,lam,dist,Tmeas,1e6*t,sum(y(:,1:2,tst),3),sum(y(:,3:end,tst),3),expflag);
    ax=axis;axis([1e6*min(t) 1e6*max(t) ax(3:4)]); xlabel('time [\mus]'); ylabel('correlation');
    s{1} = ['\itR\rm = (' mnum2str(mean(rad(1,:)),3,1) '\pm' mnum2str(std(rad(1,:)),3,1) ') nm'];
    if ~isempty(triplet)
        s{end+1} = '';
        if size(triplet,1)==1
            s{end+1} = ['\tau = ' int2str(round(mean(triplet))) ' \mus'];
        else
            tmp = ['\tau = (' int2str(round(mean(triplet(1,:))))];
            for k=2:size(triplet,1)
                tmp = [tmp ',' int2str(round(mean(triplet(k,:))))];
            end
            tmp = [tmp  ') \mus'];
            s{end+1} = tmp;
        end
    end
    s{end+1} = '';
    s{end+1} = ['\itw\rm_0 = ' int2str(round(mean(w0))) ' nm'];
    s{end+1} = ['\ita\rm_0 = ' int2str(round(mean(a0))) ' nm'];
    text(0.6,0.95,s,'VerticalAlignment','Top','units','normal')
    subplot(4,1,4)
    semilogx(1e6*t,sum(y(:,:,tst),3)-z)
    ax=axis;axis([1e6*min(t) 1e6*max(t) -max(abs(ax(3:4))) max(abs(ax(3:4)))]); xlabel('time [\mus]'); ylabel('residuals');
    set(gcf,'Position',[20 50 900 850],'PaperPositionMode','auto')
else
    y = sum(y,3);
    p = Simplex('SizeFCS',p,pmin,pmax,[],[],av,lam,dist,Tmeas,1e6*t,y(:,1:2),y(:,3:end),expflag);
    subplot(4,1,1:3)
    [err, c, z] = SizeFCS(p,av,lam,dist,Tmeas,1e6*t,y(:,1:2),y(:,3:end),expflag);
    if length(p)==3+expflag
        rad = p(1);
        w0 = p(2);
        a0 = p(3);
    else
        rad = p(1:2);
        w0 = p(3);
        a0 = p(4);
    end
    ax=axis;axis([1e6*min(t) 1e6*max(t) ax(3:4)]); xlabel('time [\mus]'); ylabel('correlation');
    if expflag>0
        triplet = p(end-expflag+1:end);
    else
        triplet = [];
    end
	s{1} = ['\itR\rm = ' mnum2str(rad(1),3,1) ' nm'];
    if ~isempty(triplet)
        s{end+1} = '';
        if length(triplet)==1
            s{end+1} = ['\tau = ' int2str(round(triplet(1))) ' \mus'];
        else
            tmp = ['\tau = (' int2str(round(triplet(1)))];
            for k=2:length(triplet)
                tmp = [tmp ',' int2str(round(triplet(k)))];
            end
            tmp = [tmp  ') \mus'];
            s{end+1} = tmp;
        end
    end
    s{end+1} = '';
    s{end+1} = ['\itw\rm_0 = ' int2str(round(w0)) ' nm'];
    s{end+1} = ['\ita\rm_0 = ' int2str(round(a0)) ' nm'];
    text(0.6,0.95,s,'VerticalAlignment','Top','units','normal')    
    subplot(4,1,4)
    semilogx(1e6*t,y-z)
    ax=axis;axis([1e6*min(t) 1e6*max(t) -max(abs(ax(3:4))) max(abs(ax(3:4)))]); xlabel('time [\mus]'); ylabel('residuals');
    pos = get(0,'ScreenSize');
    set(gcf,'Position',[0.1*pos(3) 0.05*pos(4) 0.65*pos(3) 0.85*pos(4)],'PaperPositionMode','auto')
end


