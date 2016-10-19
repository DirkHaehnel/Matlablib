% program Atto655Sat

path = 'D:\Joerg\Doc\Fcs\2Focus\Saturation\2006-03-14\';
names = dir([path 'Atto*W.mat']);
lam = [640 670]/1.333; 
av = 100e3/60;

if 1
    global pd

    for j=1:length(names)
        pow(j) = str2num(names(j).name(9:12));
    end
    for j=1:length(names)
        [y,t] = FCSCrossRead([path names(j).name]);
        eval(['load ' path names(j).name ' parameters']);
        pd = 1e-14/1e-6./Diff2Rad(6,str2num(parameters.Temperature));
%         p1(:,j) = simplex('Gauss2Fcs',[400 400 100],[],[],[],[],av,lam,[],1e6*t,y(:,1),[],[],[pd;pd]);
%         p2(:,j) = simplex('Gauss2Fcs',[400 400 100],[],[],[],[],av,lam,[],1e6*t,y(:,2),[],[],[pd;pd]);        
        p1a(:,j) = simplex('GaussFcs',[400 100],[],[],[],[],av,lam,[],1e6*t,y(:,1),[],[],[pd;pd]);
        p2a(:,j) = simplex('GaussFcs',[400 100],[],[],[],[],av,lam,[],1e6*t,y(:,2),[],[],[pd;pd]);        
    end

    save Atto655SatDiameter names pow p1a p2a -append
end

if 0
    load Atto655SatDiameter
    [pow, ind] = sort(pow);
    pow = pow/10;
    p1 = p1(:,ind);
    p2 = p2(:,ind);
    pow(1)=[]; p1(:,1)=[]; p2(:,1)=[];
    w0 = mean(p(2,:)); a0 = mean(p(3,:));
    pow0 = [10 20 30 40 50 60 70 80 90 160 320 640];
    pp = DertiOptSat(w0,a0,pow0);
    pp1 = DertiOptSat(w0,a0,pow0,2);
    save Atto655SatDiameter pp pp1 w0 a0 -append
end


