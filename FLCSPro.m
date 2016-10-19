% program FLCSPro

if 0
    path = 'D:\Daten\Gregor\050823AttoAlexaFLCS\';
    alexa = 'Alexa647_2_';
    atto = 'Atto655_';
    mix = 'Mix_';

    [resAlexa, headAlexa] = pt3Proj([path alexa '1.pt3'], [path alexa '2.pt3'], 'xfcs');
    [resAtto, headAtto] = pt3Proj([path atto '1.pt3'], [path atto '2.pt3'], 'xfcs');

    tau = resAlexa.tau';
    ind = tau>2 & tau<15;

    [resMix, headMix] = pt3Proj([path mix '1.pt3'], [path mix '2.pt3'], 'flcs', [tau(ind) tau(ind) resAlexa.tcspc(ind,:) resAtto.tcspc(ind,:) ind(ind) ind(ind)]);
end

if 1
    % load D:\Joerg\Doc\Sykora\2GFPl.mat
    % D:\MATLAB\GFPl_filter.mat

    name = {'R:\051208\hGFPtril_01.pt3','R:\051208\hGFPtril_02.pt3'};

    res = pt3pro(name{1},name{2},'tcspc');
    ind = res.tau>=3.5 & res.tau<=24;
    tau1 = simplex('expfun',[2 5],[],[],[],[],res.tau(ind),res.tcspc(ind,1));
    tau2 = simplex('expfun',[2 5],[],[],[],[],res.tau(ind),res.tcspc(ind,2));
    [err, ctau1, zz1] = expfun(tau1,res.tau(ind),res.tcspc(ind,1));
    [err, ctau2, zz2] = expfun(tau2,res.tau(ind),res.tcspc(ind,2));
    param = [res.tau(ind)', zz1*ctau1(2:3) zz2*ctau2(2:3) ones(sum(ind),2)];
    weight = res.tcspc(ind,:);
    flcs = pt3pro(name{1},name{2},'flcs',param,weight);
end

