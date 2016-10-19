% program Cy5AttoFCS

path = 'D:\Daten\Dertinger\Cy5-Atto655\';

if 0
    names = dir([path '*.t3r'])
    for j=1:length(names)
        s = names(j).name;
        [res, head] = TwoFocus2FCS([path s]);
        t1 = findstr(s,'_');
        t2 = findstr(s,'C'); t2=t2(end);
        t1 = t1(t1<t2 & t1([2:end end])>t2);
        head.Temperature = str2num(s(t1+1:t2-1));
        t1 = findstr(s,'_');
        t2 = findstr(s,'uW'); t2=t2(end);
        t1 = t1(t1<t2 & t1([2:end end])>t2);
        head.Power = str2num(s(t1+1:t2-1));
        eval(['save ' path s(1:end-4) '_joerg res head'])
    end
end

if 0
    namesatto = dir([path 'Atto655_*_parts.mat']);
    powatto = [5 5 5 10 20 30 40 50 60];
    global pd
    for j=1:length(namesatto)
        pd = 0.02;
        [dcatto(:,j) vatto(:,j) catto(:,j) w0atto(:,j) a0atto(:,j)] = FCSFit([path namesatto(j).name],[400 150],[],1);
    end
    eval(['save ' path 'AttoCy5FCS namesatto dcatto powatto vatto catto w0atto a0atto'])
end

if 0
    powcy5 = [5 10 20 30 40];
    namescy52f = dir([path '2F_Cy5_*_parts.mat']);
    global pd
    for j=1:length(namescy52f)
        pd = [0.002 10 10 10];
        [dccy52f(:,j) vcy52f(:,j) ccy52f(:,j) w0cy52f(:,j) a0cy52f(:,j), tcy52f(:,:,j)] = FCSFit([path namescy52f(j).name],[400 150],[],1);
    end
    eval(['save ' path 'AttoCy5FCS powcy5 namescy52f dccy52f vcy52f ccy52f w0cy52f a0cy52f tcy52f -append']);
end

if 0
    namescy52f = dir([path '2F_Cy5_*_parts.mat']);
    global pd
    for j=1:length(namescy52f)
        pd = [0.002 10];
        [dccy5(:,j) vcy5(:,j) ccy5(:,j) w0cy5(:,j) a0cy5(:,j), tcy5(:,j)] = FCSFit([path namescy52f(j).name],[400 150],[],1);
    end
    eval(['save ' path 'AttoCy5FCS dccy5 vcy5 ccy5 w0cy5 a0cy5 tcy5 -append']);
end

if 0
    namescy51f = dir([path '1F_Cy5_*_parts.mat']);
    global pd
    av = 100e3/58;
    lam = [640 670]/1.33;
    for j=1:length(namescy51f)
        pd = [0.002 10];
        [y t] = FCSCrossRead([path namescy51f(j).name]);
        y = sum(sum(y,3),2);
        para = simplex('GaussFcs',[400 150],[400 0],[400 inf],[],[],av,lam,[],1e6*t,y,[],1);
        a0cy51f(j) = para(2);
        dccy51f(j) = 1e-14/1e-6./pd(1);
        tcy51f(j) = pd(2);
    end
    eval(['save ' path 'AttoCy5FCS namescy51f dccy51f a0cy51f tcy51f -append']);
end

if 0
    namesatto = dir([path 'Atto655_*_joerg.mat']);
    powatto = [5 5 5 10 20 30 40 50 60];
    global pd
    for j=1:length(namesatto)
        pd = 0.02;
        [dcatto{j} vatto{j} catto{j} w0atto{j} a0atto{j}] = FCSFit([path namesatto(j).name],[400 150],[],1);
    end
    eval(['save ' path 'AttoCy5FCS_joerg namesatto dcatto powatto vatto catto w0atto a0atto'])
end

if 0
    powcy5 = [5 10 20 30 40 50 60];
    namescy52f = dir([path '2F_Cy5_*_joerg.mat']);
    global pd
    for j=1:length(namescy52f)
        pd = [0.002 10 10 10];
        [dccy52f{j} vcy52f{j} ccy52f{j} w0cy52f{j} a0cy52f{j}, tcy52f{j}] = FCSFit([path namescy52f(j).name],[400 150],[],1);
    end
    eval(['save ' path 'AttoCy5FCS_joerg powcy5 namescy52f dccy52f vcy52f ccy52f w0cy52f a0cy52f tcy52f -append']);
end

if 0
    namescy52f = dir([path '2F_Cy5_*_joerg.mat']);
    global pd
    for j=1:length(namescy52f)
        pd = [0.002 10];
        [dccy5{j} vcy5{j} ccy5{j} w0cy5{j} a0cy5{j}, tcy5{j}] = FCSFit([path namescy52f(j).name],[400 150],[],1);
    end
    eval(['save ' path 'AttoCy5FCS_joerg dccy5 vcy5 ccy5 w0cy5 a0cy5 tcy5 -append']);
end

if 0
    namescy51f = dir([path '1F_Cy5_*_joerg.mat']);
    global pd
    av = 100e3/58;
    lam = [640 670]/1.33;
    for j=1:length(namescy51f)
        pd = [0.002 10];
        [y t] = FCSCrossRead([path namescy51f(j).name]);
        y = sum(sum(y,3),2);
        para = simplex('GaussFcs',[400 150],[400 0],[400 inf],[],[],av,lam,[],1e6*t,y,[],1);
        a0cy51f(j) = para(2);
        dccy51f(j) = 1e-14/1e-6./pd(1);
        tcy51f(j) = pd(2);
    end
    eval(['save ' path 'AttoCy5FCS_joerg namescy51f dccy51f a0cy51f tcy51f -append']);
end

if 0
    clear all
    path = 'D:\Joerg\Doc\Fcs\2Focus\Saturation\2006-03-14\'
    namesatto = dir([path 'Atto655*W.mat']);
	for j=1:length(namesatto)
       pow(j) = str2num(namesatto(j).name(end-9:end-6))/10;
    end
    av = 100e3/58;
    lam = [640 670]/1.33;
    distv = [380 390 400];
    global pd
    for k=1:length(distv)
        para = [av lam distv(k)];
        for j=1:length(namesatto)
            pd = 0.02;
            [dcatto(j,k) vatto(j,k) catto(j,k) w0atto(j,k) a0atto(j,k) tmp c(:,:,j,k) err(j,k)] = FCSFit([path namesatto(j).name],[300 100],[],[],para);
            eval(['print -dpng ' namesatto(j).name(1:end-4) '_' mint2str(distv(k),3)])
        end
        eval(['save ' path 'AttoFCS_dist namesatto dcatto vatto catto w0atto a0atto c err'])
    end
end

if 1
    clear all
    path = 'D:\Joerg\Doc\Fcs\2Focus\Saturation\2006-03-14\'
    load D:\Joerg\Doc\Fcs\2Focus\Saturation\2006-03-14\AttoFCS_dist.mat dcatto
    namesatto = dir([path 'Atto655*W.mat']);
	for j=1:length(namesatto)
       pow(j) = str2num(namesatto(j).name(end-9:end-6))/10;
    end
    av = 100e3/58;
    lam = [640 670]/1.33;
    distv = [380 390 400];
    global pd
    for k=2:length(distv)
        para = [av lam distv(k)];
        for j=1:length(namesatto)
            pd = 1e-8/dcatto(j,k);
            [y t] = FCSCrossRead([path namesatto(j).name]);
            y = sum(sum(y,3),2);
            para = simplex('GaussFcs',[400 150],[0 0],[inf inf],[],[],av,lam,[],1e6*t,y,[],[],[pd;pd]);
            w0atto1f(j,k) = para(1);
            a0atto1f(j,k) = para(2);
        end
        eval(['save ' path 'AttoFCS_mono pow namesatto dcatto w0atto1f a0atto1f'])
    end
end
