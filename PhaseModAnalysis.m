% PhaseModAnalysis

names = {'f:\Joerg\Doc\Fcs\2Focus\PhaseModulation\AT655 20uW Sin_10us.spc', ...
    'c:\Joerg\Doc\Fcs\2Focus\PhaseModulation\AT655 80uW Sin_10us again.spc', ...
    'c:\Joerg\Doc\Fcs\2Focus\PhaseModulation\AT655 80uW Sin_10us.spc', ...
    'c:\Joerg\Doc\Fcs\2Focus\PhaseModulation\AT655 200uW Sin_10us.spc', ...
    'c:\Joerg\Doc\Fcs\2Focus\PhaseModulation\AT655 160uW Sin_10us.spc'};

alv = 3:0.1:4;
for j=1:1%length(names)
    name = names{j};
    [res, head] = PhaseMod2FCS(name, 10e-6);

    global pd
    pd = 0.1;

    ind = res.autotime<1e-1;
    data.t = res.autotime(ind);
    for k=1:length(alv)
        al = alv(k);
        a = squeeze(sum(res.auto1(ind,1,2,:)+res.auto1(ind,2,1,:)+al*res.auto2(ind,1,2,:)+al*res.auto2(ind,2,1,:),4));
        b = squeeze(sum(res.auto1(ind,1,2,:)+res.auto1(ind,2,1,:)-al*res.auto2(ind,1,2,:)-al*res.auto2(ind,2,1,:),4));
        data.y = permute(cat(3, a, a, b, b),[1 3 2]);
        [dc(k) bla bla w0(k) a0(k) bla bla bla err(k)] = FCSFit(data,[350, 150],[],1,[100e3/60 [635 670]/1.33 433]);
    end

    eval(['print -dpng -r300 ''' name(1:end-4) ''''])
end

return

al = 3.5;
a=squeeze(sum(res.auto1(:,1,2,:)+res.auto1(:,2,1,:)+al*res.auto2(:,1,2,:)+al*res.auto2(:,2,1,:),4));
b=squeeze(sum(res.auto1(:,1,2,:)+res.auto1(:,2,1,:)-al*res.auto2(:,1,2,:)-al*res.auto2(:,2,1,:),4));
data.t = res.autotime;
data.y = permute(cat(3, a, a, b, b),[1 3 2]);
[dc v conc w0 a0] = FCSFit(data,[350, 150],[],1,[100e3/60 [635 670]/1.33 433]);