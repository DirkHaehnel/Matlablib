function [dc w0 a0 triplet c velo err z] = FCSFit_MHC(name,p,expflag,bootstrap,para,bounds,pmin,pmax)

% Function FCSFit for fitting 2fFCS data
% (c) Joerg Enderlein, http://www.joerg-enderlein.de (2010)

close all
global pd

if nargin<2 || isempty(p)
    p = [400 150];
end

if nargin>4 && ~isempty(para)
    av = para(1);
    lam = para(2:3);
    dist = para(4);
else
    av = 75e3/60;
    lam = [570 670]/1.33;
    dist = 394;
end

if length(p)==3 && p(3)==1i
    p = p(1:2);
    flag2d = 1;
else
    flag2d = 0;
    if length(p)>2
        p(3:end) = pd(1)*p(3:end)/1e3;
    end
end

if isstruct(name)
    if isfield(name,'auto')
        [y,t] = FCSCrossRead(name);
    else
        y = name.y;
        t = name.t;
    end
else
    [y,t] = FCSCrossRead(name);
end

if size(y,3)>1
    y = y(:,:,~isnan(sum(sum(y,1),2)));
end

if isempty(pd)
    disp('Please define global variable pd');
    return
end
if nargin<3 || isempty(expflag)
    expflag = [length(pd)-1 0];
elseif length(expflag)==1
    expflag = [expflag 0];
end
if nargin<6 || isempty(bounds)
    bounds = [];
end
if nargin<7 || isempty(pmin)
    pmin = 0*p;
    if length(p)>2
        pmin(3:end) = -inf;
    end
end
if nargin<8 || isempty(pmax)
    pmax = inf*ones(size(p));
end
if nargin>3 && ~(isempty(bootstrap) || (bootstrap==0)) && size(y,3)>1
    if size(y,3)>1
        weight = 1./std(y,0,3);
    else
        weight = [];
    end
    p0 = Simplex('GaussFCS',p,pmin,pmax,[],[],av,lam,dist,1e6*t,mean(y(:,1:2,:),3),mean(y(:,3:end,:),3),expflag,bounds,flag2d,weight);
    pd0 = pd;
    for j=1:bootstrap
        pd = pd0;
        ind = randi(size(y,3),[1 size(y,3)]);
        p = Simplex('GaussFCS',p0,pmin,pmax,[],[],av,lam,dist,1e6*t,mean(y(:,1:2,ind),3),mean(y(:,3:end,ind),3),expflag,bounds,flag2d,weight);
        [err(j), c(:,:,j)] = GaussFCS(p,av,lam,dist,1e6*t,mean(y(:,1:2,ind),3),mean(y(:,3:end,ind),3),expflag,bounds,flag2d,weight);
        if expflag(1)>0
            triplet(:,j) = pd(end-expflag+1:end);
        else 
            triplet = [];
        end
        if length(p)>2
            velo(:,j) = zeros(3,1);
            velo(1:length(p)-2,j) = p(3:end)/pd(1)*1e3;
            p(3:end) = [];
        else
            velo = [];
        end
        if flag2d
            w0(1,j) = p(1);
            w0(2,j) = p(2);
            p(2) = sqrt(sum(p(1:2).^2)/2); p(1) = [];
        else
            w0(j) = p(1);
        end
        if length(p)>1
            a0(j) = p(2);
        else
            a0(j) = 0;
        end
        dc(1:length(pd)-sum(expflag),j) = 1e-14/1e-6./pd(1:end-sum(expflag));
        if expflag(2)>0
            dc(length(pd)-sum(expflag)+1:length(pd)-sum(expflag)+expflag(2),j) = pd(end-sum(expflag)+1:end-expflag(1));
        end
    end
    subplot(4,1,1:3)
    [tmp, tmp, z] = GaussFCS(p0,av,lam,dist,1e6*t,mean(y(:,1:2,:),3),mean(y(:,3:end,:),3),expflag,bounds,flag2d);
    ax=axis;axis([1e6*min(t) 1e6*max(t) ax(3:4)]); ylabel('correlation / (\ith\rm\nu)^2/s^2'); 
    for j=1:length(pd)-sum(expflag)
        tmp = -floor(log10(mean(dc(j,:))));
        s{j} = ['\itD\rm_' mint2str(j) ' = (' mnum2str(mean(dc(j,:))*10^tmp,1,2)  '\pm' mnum2str(std(dc(j,:))*10^tmp,1,2) ')\cdot10^{-' int2str(tmp) '} cm^2/s'];
    end
    if ~isempty(velo)
        if ~(mean(velo(1,:))==0)
            s{end+1} = ['\itv_x\rm = (' mnum2str(mean(velo(1,:)),1,2) '\pm' mnum2str(std(velo(1,:)),1,2) ') \mum/s'];
        end
        if ~(mean(velo(2,:))==0)
            s{end+1} = ['\itv_y\rm = (' mnum2str(mean(velo(2,:)),1,2) '\pm' mnum2str(std(velo(2,:)),1,2) ') \mum/s'];
        end
        if ~(mean(velo(3,:))==0)
            s{end+1} = ['\itv_z\rm = (' mnum2str(mean(velo(3,:)),1,2) '\pm' mnum2str(std(velo(3,:)),1,2) ') \mum/s'];
        end
    end
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
    if a0>0
        s{end+1} = ['\ita\rm_0 = ' int2str(round(mean(a0))) ' nm'];
    end
    text(0.6,0.95,s,'VerticalAlignment','Top','units','normal')
    subplot(4,1,4)
    semilogx(1e6*t,mean(y,3)-z)
    ax=axis;axis([1e6*min(t) 1e6*max(t) -max(abs(ax(3:4))) max(abs(ax(3:4)))]); xlabel('time [\mus]'); ylabel('residuals');
    set(gcf,'Position',[20 50 900 850],'PaperPositionMode','auto')
else
    if size(y,3)>1
        weight = 1./std(y,0,3);
    else
        weight = [];
    end
    y = mean(y,3);
    p = Simplex('GaussFcs',p,pmin,pmax,[],[],av,lam,dist,1e6*t,y(:,1:2),y(:,3:end),expflag,bounds,flag2d,weight);
    subplot(4,1,1:3)
    [err, c, z] = GaussFCS(p,av,lam,dist,1e6*t,y(:,1:2),y(:,3:end),expflag,bounds,flag2d,weight);
    if length(p)>2
        velo = zeros(3,1);
        velo(1:length(p)-2) = p(3:end)/pd(1)*1e3;
        p(3:end) = [];
    else
        velo = [];
    end
    if flag2d
        w0(1) = p(1);
        w0(2) = p(2);
        p(2) = sqrt(sum(p(1:2).^2)/2); p(1) = [];
    else
        w0 = p(1);
    end
    if length(p)>1
        a0 = p(2);
    else
        a0 = 0;
    end
    dc = 1e-14/1e-6./pd(1:end-sum(expflag));
    if expflag(2)>0
        dc(end+1:end+expflag(2)) = pd(end-sum(expflag)+1:end-expflag(1));
    end
    ax=axis;axis([1e6*min(t) 1e6*max(t) ax(3:4)]); ylabel('correlation / (\ith\rm\nu)^2/s^2'); 
    if expflag(1)>0
        triplet = pd(end-expflag+1:end);
    else
        triplet = [];
    end
    for j=1:length(pd)-sum(expflag)
        tmp = -floor(log10(dc(j)));
        s{j} = ['\itD\rm_' mint2str(j) ' = ' mnum2str(dc(j)*10^tmp,1,2) '\cdot10^{-' int2str(tmp) '} cm^2/s'];
    end
    if ~isempty(velo)
        s{end+1} = ['\bfv\rm = (' mnum2str(velo(1),1,2) ', ' mnum2str(velo(2),1,2) ', ' mnum2str(velo(3),1,2) ') \mum/s'];
    end
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
    if a0>0
        s{end+1} = ['\ita\rm_0 = ' int2str(round(a0)) ' nm'];
    end
    text(0.6,0.95,s,'VerticalAlignment','Top','units','normal')    
    subplot(4,1,4)
    semilogx(1e6*t,y-z)
    ax=axis;axis([1e6*min(t) 1e6*max(t) -max(abs(ax(3:4))) max(abs(ax(3:4)))]); xlabel('time [\mus]'); ylabel('residuals');
    pos = get(0,'ScreenSize');
    set(gcf,'Position',[0.1*pos(3) 0.05*pos(4) 0.65*pos(3) 0.85*pos(4)],'PaperPositionMode','auto')
end

