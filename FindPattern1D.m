function [err, bdata, cdata, sdata, xc, bc, cc, sc, len] = FindPattern1D(data,mask,sup,bck,sze,tsh,fun,flag)

if nargin<7 || isempty(fun)
    fun = 'cdata>tsh*sqrt(err)';
end
n1 = size(mask,1);
n2 = size(mask,2);
if nargin<3 || isempty(bck)
    bck = ones(n1,n2);
elseif size(bck,2)<size(mask,2)
    bck = repmat(bck,[1 n2]);
end
if nargin<4 || isempty(sup)
    sup = ones(n1,n2);
elseif size(sup,2)<size(mask,2)
    sup = repmat(sup,[1 n2]);
end
if nargin<8 || isempty(flag)
    flag = true;
end
for j=1:n2
    bck(:,j) = bck(:,j)/sqrt(sum(sup(:,j).*bck(:,j).^2));
    mask(:,j) = mask(:,j).*(sup(:,j)>0);
    mask(:,j) = mask(:,j)./sqrt(sum(sup(:,j).*mask(:,j).^2));
end
err = inf*data;
cdata = data;
bdata = data;
sdata = 1+0*data;
for s = 1:n2
    data0 = mConv(data,sup(:,s).*bck(:,s));
    data02 = mConv(data.^2,sup(:,s));
    crs = sum(sup(:,s).*bck(:,s).*mask(:,s));
    crs = inv([1 crs; crs 1]);
    data1 = mConv(data,sup(:,s).*mask(:,s));
    tmperr = data02 - crs(1,1)*data0.*data0 - 2*crs(1,2)*data0.*data1 - crs(2,2)*data1.*data1;
    tmpbdata = crs(1,1)*data0 + crs(1,2)*data1;
    tmpcdata = crs(1,2)*data0 + crs(2,2)*data1;
    cdata(tmperr<err) = tmpcdata(tmperr<err);
    bdata(tmperr<err) = tmpbdata(tmperr<err);
    sdata(tmperr<err) = s;
    err(tmperr<err) = tmperr(tmperr<err);
end

if nargout>4
    if nargin<5 || isempty(sze)
        sze = 1;
    end
    if nargin<6 || isempty(tsh)
        tsh = 1;
    end
    eval(['[cl,len]=mCluster(' fun ',sze);'])
    a = 1:length(data);
    cnt = 1; xc = []; bc = []; cc = []; sc = [];
    for j=1:length(len)
        ind = cl==j;
        tmpxc = a(ind);
        tmperr = err(ind);
        tmpb = bdata(ind);
        tmpc = cdata(ind);
        tmps = sdata(ind);
        ind = tmperr==min(tmperr);
        sc(cnt) = tmps(ind);
        bc(cnt) = tmpb(ind);
        cc(cnt) = tmpc(ind);
        xc(cnt) = tmpxc(ind);
        cnt = cnt+1;
    end
    if ~isempty(xc)
        xc = round(xc);
        if flag
            ind = xc<(n1+1)/2 | xc>length(data)-(n1-1)/2;
            xc(ind)=[]; bc(ind) = []; cc(ind) = []; sc(ind) = []; len(ind) = [];
        end
    end
end
