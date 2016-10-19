function [err, bim, cim, sim, xc, yc, bc, cc, sc, len, cl] = FindOrientation(im, r, p, sze, tsh, fun, f00, f01, f02, f10, f11, f12, f20, f21, f22, flag);

if nargin<6 | isempty(fun)
    fun = 'cim>tsh*sqrt(err)';
end
n1 = size(f00,1);
n2 = size(f00,2);
if nargin<16 | isempty(flag)
    flag = true;
end

err = inf*im;
cim = im;
bim = im;
sim = 1+0*im;
bck = disk(round((size(f00)-1)/2));
imb = mconv2(im,bck);
im00 = mconv2(im,f00);
im01 = mconv2(im,f01);
im02 = mconv2(im,f02);
im10 = mconv2(im,f10);
im11 = mconv2(im,f11);
im12 = mconv2(im,f12);
im20 = mconv2(im,f20);
im21 = mconv2(im,f21);
im22 = mconv2(im,f22);

tmp = sum(sum(bck.*im00));
cr00 = inv([1 tmp; tmp 1]);
tmp = sum(sum(bck.*im01));
cr01 = inv([1 tmp; tmp 1]);
tmp = sum(sum(bck.*im02));
cr02 = inv([1 tmp; tmp 1]);
tmp = sum(sum(bck.*im10));
cr10 = inv([1 tmp; tmp 1]);
tmp = sum(sum(bck.*im11));
cr11 = inv([1 tmp; tmp 1]);
tmp = sum(sum(bck.*im12));
cr12 = inv([1 tmp; tmp 1]);
tmp = sum(sum(bck.*im20));
cr20 = inv([1 tmp; tmp 1]);
tmp = sum(sum(bck.*im21));
cr21 = inv([1 tmp; tmp 1]);
tmp = sum(sum(bck.*im22));
cr22 = inv([1 tmp; tmp 1]);

    im1 = mconv2(im,squeeze(mask(:,:,s)));
    tmperr = im02 - crs(1,1)*im0.*im0 - 2*crs(1,2)*im0.*im1 - crs(2,2)*im1.*im1;
    tmpcim = crs(1,2)*im0 + crs(2,2)*im1;
    tmpbim = crs(1,1)*im0 + crs(1,2)*im1;
    cim(tmperr<err) = tmpcim(tmperr<err);
    bim(tmperr<err) = tmpbim(tmperr<err);
    sim(tmperr<err) = s;
    err(tmperr<err) = tmperr(tmperr<err);
    s
    end
end

    
    
if nargout>4
    if nargin<4 | isempty(sze)
        sze = 1;
    end
    if nargin<5 | isempty(tsh)
        tsh = 1;
    end
    eval(['[cl,len]=mCluster(' fun ',sze);'])
    [a,b] = meshgrid(1:size(im,2),1:size(im,1));
    cnt = 1; xc = []; yc = []; bc = []; cc = []; sc = [];
    for j=1:length(len)
        ind = cl==j;
        tmpxc = a(ind);
        tmpyc = b(ind);
        tmperr = err(ind);
        tmpb = bim(ind);
        tmpc = cim(ind);
        tmps = sim(ind);
        ind = tmperr==min(tmperr);
        sc(cnt) = tmps(ind);
        bc(cnt) = tmpb(ind);
        cc(cnt) = tmpc(ind);
        xc(cnt) = tmpxc(ind);
        yc(cnt) = tmpyc(ind);    
        cnt = cnt+1;
    end
    if ~isempty(xc)
        xc = round(xc);
        yc = round(yc);
        if flag
            ind = xc<(n2+1)/2 | xc>size(im,2)-(n2-1)/2 | yc<(n1+1)/2 | yc>size(im,1)-(n1-1)/2; 
            xc(ind)=[]; yc(ind)=[]; bc(ind) = []; cc(ind) = []; sc(ind) = []; len(ind) = [];
        end
    end
end

