function [rout, nout, pxout, pyout, Lambda] = Refraction(rin, nin, pxin, pyin, R, z, n)

if R==0
   Lambda = (z - rin(:,3))./nin(:,3);
else
   Lambda = -sum(rin.*nin,2) + nin(:,3)*z;
   Lambda = Lambda*[1 1] + (sign(R)*sqrt(Lambda.^2 + R^2 - sum(rin.^2,2) + 2*rin(:,3)*z - z^2))*[1 -1]; 
   Lambda(Lambda<0) = inf;
   Lambda = min(Lambda,[],2);
end
rout = rin + (Lambda*[1 1 1]).*nin;

if R==0
   direx = ones(size(rin,1),1)*[1 1 1];
else
   direx = [rout(:,1) rout(:,2) rout(:,3)-z];
   direx = direx./(sqrt(sum(direx.^2,2))*[1 1 1]);
end
tmp = sum(direx.*nin,2);
nout = 1/n*(nin - (tmp*[1 1 1]).*direx) + direx.*((sign(tmp).*sqrt(1-1/n^2*(1-sum(direx.*nin,2).^2)))*[1 1 1]);

if nargout>2 && ~isempty(pxin)
    ts = cross(nout,nin);
    ww = sqrt(sum(ts.^2,2));
    tst = ww>1e-5;
    ts(tst,:) = ts(tst,:)./(ww(tst)*[1 1 1]);
    tp = cross(ts(tst,:),nin(tst,:));
    tpprime = cross(ts(tst,:),nout(tst,:));

    pxout = pxin;
    pxout(tst,:) = (sum(pxin(tst,:).*ts(tst,:),2)*[1 1 1]).*ts(tst,:) + (sum(pxin(tst,:).*tp,2)*[1 1 1]).*tpprime;
    if ~isempty(pyin)
        pyout = pyin;
        pyout(tst,:) = (sum(pyin(tst,:).*ts(tst,:),2)*[1 1 1]).*ts(tst,:) + (sum(pyin(tst,:).*tp,2)*[1 1 1]).*tpprime;
    else
        pyout = [];
    end
end