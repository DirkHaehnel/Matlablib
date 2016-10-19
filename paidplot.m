function paidplot(res, x, y);

if nargin==1
    if size(res.paid,3)>1
        pcolor(res.autotime,res.nums,log10(sum(res.paid,3)'));
    else
        pcolor(res.autotime,res.nums,log10(res.paid'));
    end
else if nargin==3
        pcolor(res,x,log10(y'));
    end
end
set(gca,'xscale','log','yscale','log')
shading flat
times16
xlabel('time [s]')
ylabel('photon number')