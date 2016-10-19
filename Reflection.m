function [rout, nout, pxout, pyout, Lambda] = Reflection(rin, nin, pxin, pyin, R, z, n)

if R==0
   Lambda = (z - rin(:,3))./nin(:,3);
else
   Lambda = -sum(rin.*nin,2) + nin(:,3)*z;
   Lambda = Lambda + sign(R)*sqrt(Lambda.^2 + R^2 - sum(rin.^2,2) + 2*rin(:,3)*z - z^2);
end
rout = rin + (Lambda*[1 1 1]).*nin;

if R==0
   direx = ones(size(rin,1),1)*[1 1 1];
else
   direx = [rout(:,1) rout(:,2) rout(:,3)-z];
   direx = direx./(sqrt(sum(direx.^2,2))*[1 1 1]);
end
w = sum(direx.*nin,2);
nout = nin - 2*(w*[1 1 1]).*direx;

if nargout>2 && ~isempty(pxin)
    [rp, rs] = Fresnel(w,1,n);

    ts = cross(nout,nin);
    ww = sqrt(sum(ts.^2,2));
    tst = ww>1e-5;
    ts(tst,:) = ts(tst,:)./(ww(tst)*[1 1 1]);
    tp = cross(ts(tst,:),nin(tst,:));
    tpprime = cross(ts(tst,:),nout(tst,:));

    pxout = pxin;
    pxout(tst,:) = ((rs.'.*sum(pxin(tst,:).*ts(tst,:),2))*[1 1 1]).*ts(tst,:) + ((rp.'.*sum(pxin(tst,:).*tp,2))*[1 1 1]).*tpprime;
    if ~isempty(pyin)
        pyout = pyin;
        pyout(tst,:) = ((rs.'.*sum(pyin(tst,:).*ts(tst,:),2))*[1 1 1]).*ts(tst,:) + ((rp.'.*sum(pyin(tst,:).*tp,2))*[1 1 1]).*tpprime;
    else
        pyout = [];
    end
end

