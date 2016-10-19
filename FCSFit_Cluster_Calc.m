function [p, err, c, triplet, velo, w0, a0, dc z]= FCSFit_Cluster_Calc(file_name)
% This program is the helper for FCSFit_Cluster.m which is called while
% each job is carried out

% (c) Narain Karedla 2013

   load (file_name)
        global pd
        t=databs.t;
        y=databs.y;
        pd=pd0;
        p = Simplex('GaussFcs',p0,pmin,pmax,[],[],av,lam,delta,1e6*t,y(:,1:2),y(:,3:end),expflag,bounds,flag2d,weight,bld);
        [err, c, z] = GaussFcs(p,av,lam,delta,1e6*t,y(:,1:2),y(:,3:end),expflag,bounds,flag2d,weight,bld);
        if expflag(1)>0
            triplet = pd(end-expflag+1:end);
        else 
            triplet = [];
        end
        if length(p)>2
            velo = zeros(3,1);
            velo(1:length(p)-2) = p(3:end)/pd(1)*1e3;
            p(3:end) = [];
        else
            velo = [];
        end
        if flag2d
            w0(1) = p(1);
            w0(2) = p(2);
            p(2) = sqrt(sum(p(1:2).^2)/2); p(1) = [];
        else
            w0 = p(1);
        end
        if length(p)>1
            a0 = p(2);
        else
            a0 = 0;
        end
        dc(1:length(pd)-sum(expflag)) = 1e-14/1e-6./pd(1:end-sum(expflag));
        if expflag(2)>0
            dc(length(pd)-sum(expflag)+1:length(pd)-sum(expflag)+expflag(2)) = pd(end-sum(expflag)+1:end-expflag(1));
        end
end