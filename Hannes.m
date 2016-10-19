% program for evaluating the antibunching measurements of Hannes Neuweiler

if 1
    %y(:,1)
    [x, flag, tcspc, head, num] = tttrRead('c:\JOERG\Doc\Neuweiler\1antiP122MR121P10.t3r');
    t = 5:10:head.NChannels;
    NCounts = head.NCounts;
    TAcq = head.TAcq;
    rate = sqrt(NCounts/TAcq*1e3/5e-6);
    y = mhist(max(tcspc)-tcspc,t);
    
    %y(:,2)
    [x, flag, tcspc, head, num] = tttrRead('c:\JOERG\Doc\Neuweiler\1antiP122MR121P100.t3r');
    tmp = mhist(max(tcspc)-tcspc,t);
    NCounts(2) = head.NCounts;
    TAcq(2) = head.TAcq;
    [x, flag, tcspc, head, num] = tttrRead('c:\JOERG\Doc\Neuweiler\2antiP122MR121P100.t3r');
    tmp = tmp + mhist(max(tcspc)-tcspc,t);
    NCounts(2) = NCounts(2) + head.NCounts;
    TAcq(2) = TAcq(2) + head.TAcq;
    [x, flag, tcspc, head, num] = tttrRead('c:\JOERG\Doc\Neuweiler\3antiP122MR121P100.t3r');
    tmp = tmp + mhist(max(tcspc)-tcspc,t);
    NCounts(2) = NCounts(2) + head.NCounts;
    TAcq(2) = TAcq(2) + head.TAcq;
    rate = [rate sqrt(NCounts(end)/TAcq(end)*1e3/5e-6)];
    y = [y tmp];
    
    %y(:,3)
    [x, flag, tcspc, head, num] = tttrRead('c:\JOERG\Doc\Neuweiler\1antiP122MMR121P10.t3r');
    NCounts = [NCounts head.NCounts];
    TAcq = [TAcq head.TAcq];
    rate = [rate sqrt(NCounts(end)/TAcq(end)*1e3/5e-6)];    
    y = [y mhist(max(tcspc)-tcspc,t)];

    %y(:,4)
    [x, flag, tcspc, head, num] = tttrRead('c:\JOERG\Doc\Neuweiler\1antiP122MMR121P100.t3r');
    NCounts = [NCounts head.NCounts];
    TAcq = [TAcq head.TAcq];
    rate = [rate sqrt(NCounts(end)/TAcq(end)*1e3/5e-6)];    
    y = [y mhist(max(tcspc)-tcspc,t)];
     
    %y(:,5)
    [x, flag, tcspc, head, num] = tttrRead('c:\JOERG\Doc\Neuweiler\1antiP1801MR121P10.t3r');
    NCounts = [NCounts head.NCounts];
    TAcq = [TAcq head.TAcq];
    rate = [rate sqrt(NCounts(end)/TAcq(end)*1e3/5e-6)];    
    y = [y mhist(max(tcspc)-tcspc,t)];
    
    %y(:,6)
    [x, flag, tcspc, head, num] = tttrRead('c:\JOERG\Doc\Neuweiler\1antiP1801MR121P100.t3r');
    NCounts = [NCounts head.NCounts];
    TAcq = [TAcq head.TAcq];
    rate = [rate sqrt(NCounts(end)/TAcq(end)*1e3/5e-6)];    
    y = [y mhist(max(tcspc)-tcspc,t)];
    
    %y(:,7)
    [x, flag, tcspc, head, num] = tttrRead('c:\JOERG\Doc\Neuweiler\1antiP1801MMR121P10.t3r');
    NCounts = [NCounts head.NCounts];
    TAcq = [TAcq head.TAcq];
    rate = [rate sqrt(NCounts(end)/TAcq(end)*1e3/5e-6)];    
    y = [y mhist(max(tcspc)-tcspc,t)];

    %y(:,8)
    [x, flag, tcspc, head, num] = tttrRead('c:\JOERG\Doc\Neuweiler\1antiP1801MMR121P100.t3r');
    NCounts = [NCounts head.NCounts];
    TAcq = [TAcq head.TAcq];
    rate = [rate sqrt(NCounts(end)/TAcq(end)*1e3/5e-6)];    
    y = [y mhist(max(tcspc)-tcspc,t)];

    
    %y(:,9)
    [x, flag, tcspc, head, num] = tttrRead('c:\JOERG\Doc\Neuweiler\1antiMR121P10.t3r');
    NCounts = [NCounts head.NCounts];
    TAcq = [TAcq head.TAcq];
    rate = [rate sqrt(NCounts(end)/TAcq(end)*1e3/5e-6)];    
    y = [y mhist(max(tcspc)-tcspc,t)];

    %y(:,10)
    [x, flag, tcspc, head, num] = tttrRead('c:\JOERG\Doc\Neuweiler\1antiMR121P100.t3r');
    NCounts = [NCounts head.NCounts];
    TAcq = [TAcq head.TAcq];
    rate = [rate sqrt(NCounts(end)/TAcq(end)*1e3/5e-6)];    
    y = [y mhist(max(tcspc)-tcspc,t)];
    
    t = head.Resolution*t;
    
    t(end-1:end) = [];
    y(end-1:end,:) = [];
    
    p10=simplex('expfun',1/1e3,0,[],[],[],t,y(:,9));
    p100=simplex('expfun',1/1e3,0,[],[],[],t,y(:,10))    
    kiso0 = (p100-p10)/9;
    kt0 = p10 - kiso0;

    p10=simplex('expfun',1/1e3,0,[],[],[],t,y(:,3));
    p100=simplex('expfun',1/1e3,0,[],[],[],t,y(:,4))    
    kiso122 = (p100-p10)/9;
    kt122 = p10 - kiso122;
    
    p10=simplex('expfun',1/1e3,0,[],[],[],t,y(:,7));
    p100=simplex('expfun',1/1e3,0,[],[],[],t,y(:,8))    
    kiso180 = (p100-p10)/9;
    kt180 = p10 - kiso180;

    keq122 = 0.31/(1-0.31);
    keq180 = 0.17/(1-0.17);
    
    p = [1/1e4 1/1e2 1/1e3];
    for j=1:size(y,2) yy(:,j)=y(:,j)/sum(y(:,j)); end
    p122 = simplex('HannesExpFit', p, zeros(1,length(p)), [], [], [], t, yy(:, [1 2]), kt122, 10, keq122);
    p180 = simplex('HannesExpFit', p, zeros(1,length(p)), [], [], [], t, yy(:, [5 6]), kt180, 10, keq180);    
    
    %figures
    close all
    [err, z122] = HannesExpFit(p122, t, y(:, [1 2])*diag([3 1]), kt122, 10, keq122);
    plot(t, 3*y(:,1), 'o', t, y(:,2), 'v', 'markersize', 5);
    hold
    plot(t, z122, 'linewidth', 2);
    hold
    ax = axis;
    axis([0 5e3 0 1.6e3]);
    times20
    xlabel('time [ns]');
    ylabel('photon distance distribution [a.u.]');
    title('P122');
    
    figure
    [err, z180] = HannesExpFit(p180, t, y(:, [5 6])*diag([3 1]), kt180, 10, keq180);
    plot(t, 3*y(:,5), 'o', t, y(:,6), 'v', 'markersize', 5); 
    hold
    plot(t, z180, 'linewidth', 2);
    hold
    ax = axis;
    axis([0 5e3 0 1.6e3]);
    times20
    xlabel('time [ns]');
    ylabel('photon distance distribution [a.u.]');
    title('P1801');

    figure
    p10=simplex('expfun',1/1e3,0,[],[],[],t,y(:,7));
    [err, c, z, z10] = expfun(p10,t,y(:,7));
    p100=simplex('expfun',1/1e3,0,[],[],[],t,y(:,8))    
    [err, c, z, z100] = expfun(p100,t,y(:,8));    
    plot(t, y(:,7)/2.5, 'o', t, y(:,8)/5, 'v', 'markersize', 5); 
    hold
    plot(t, z10/2.5, t, z100/5, 'linewidth', 2);
    hold
    ax = axis;
    axis([0 5e3 0 2.5e3]);
    times20
    xlabel('time [ns]');
    ylabel('photon distance distribution [a.u.]');
    title('P1801Mutant');

    figure
    p10=simplex('expfun',1/1e3,0,[],[],[],t,y(:,3));
    [err, c, z, z10] = expfun(p10,t,y(:,3));
    p100=simplex('expfun',1/1e3,0,[],[],[],t,y(:,4))    
    [err, c, z, z100] = expfun(p100,t,y(:,4));    
    plot(t, y(:,3), 'o', t, y(:,4), 'v', 'markersize', 5); 
    hold
    plot(t, z10, t, z100, 'linewidth', 2);
    hold
    ax = axis;
    axis([0 5e3 0 2.5e3]);
    times20
    xlabel('time [ns]');
    ylabel('photon distance distribution [a.u.]');
    title('P122Mutant');

    figure
    p10=simplex('expfun',1/1e3,0,[],[],[],t,y(:,9));
    [err, c, z, z10] = expfun(p10,t,y(:,3));
    p100=simplex('expfun',1/1e3,0,[],[],[],t,y(:,10))    
    [err, c, z, z100] = expfun(p100,t,y(:,4));    
    plot(t, y(:,3), 'o', t, y(:,4), 'v', 'markersize', 5); 
    hold
    plot(t, z10, t, z100, 'linewidth', 2);
    hold
    ax = axis;
    axis([0 5e3 0 2.5e3]);
    times20
    xlabel('time [ns]');
    ylabel('photon distance distribution [a.u.]');
    title('MR121');
end

if 0
    %y(:,:,1)
    [res, head] = tttrpro('c:\JOERG\Doc\Neuweiler\1P122MR121P10.t3r','fimda');
    y = res.fida;
    nums = res.nums;
    autotime = res.autotime;
    NCounts = head.NCounts;
    TAcq = head.TAcq;
    
    %y(:,2)
    [res, head] = tttrpro('c:\JOERG\Doc\Neuweiler\1P122MR121P100.t3r','fimda');
    y(:,:,2) = res.fida;
    NCounts = [NCounts head.NCounts];
    TAcq = [TAcq head.TAcq];
    
    %y(:,3)
    [res, head] = tttrpro('c:\JOERG\Doc\Neuweiler\1P122MMR121P10.t3r','fimda');
    y(:,:,3) = res.fida;
    NCounts = [NCounts head.NCounts];
    TAcq = [TAcq head.TAcq];

    %y(:,4)
    [res, head] = tttrpro('c:\JOERG\Doc\Neuweiler\1P122MMR121P100.t3r','fimda');
    y(:,:,4) = res.fida;
    NCounts = [NCounts head.NCounts];
    TAcq = [TAcq head.TAcq];
     
    %y(:,5)
    [res, head] = tttrpro('c:\JOERG\Doc\Neuweiler\1P1801MR121P10.t3r','fimda');
    y(:,:,5) = res.fida;
    NCounts = [NCounts head.NCounts];
    TAcq = [TAcq head.TAcq];
    
    %y(:,6)
    [res, head] = tttrpro('c:\JOERG\Doc\Neuweiler\1P1801MR121P100.t3r','fimda');
    y(:,:,6) = res.fida;
    NCounts = [NCounts head.NCounts];
    TAcq = [TAcq head.TAcq];
    
    %y(:,7)
    [res, head] = tttrpro('c:\JOERG\Doc\Neuweiler\1P1801MMR121P10.t3r','fimda');
    y(:,:,7) = res.fida;
    NCounts = [NCounts head.NCounts];
    TAcq = [TAcq head.TAcq];
    
    %y(:,8)
    [res, head] = tttrpro('c:\JOERG\Doc\Neuweiler\1P1801MMR121P100.t3r','fimda');
    y(:,:,8) = res.fida;
    NCounts = [NCounts head.NCounts];
    TAcq = [TAcq head.TAcq];

    %y(:,9)
    [res, head] = tttrpro('c:\JOERG\Doc\Neuweiler\1MR121P10.t3r','fimda');
    y(:,:,9) = res.fida;
    NCounts = [NCounts head.NCounts];
    TAcq = [TAcq head.TAcq];
    
    %y(:,10)
    [res, head] = tttrpro('c:\JOERG\Doc\Neuweiler\1MR121P100.t3r','fimda');
    y(:,:,10) = res.fida;
    NCounts = [NCounts head.NCounts];
    TAcq = [TAcq head.TAcq];
end