cd E:\Daten\Gregor\041216

% for j=1:length(names) [data(j) hdata(j)]=tttrpro(names(j).name,'xfcsfilter'); end
% save data0 data hdata names
% for j=1:length(names) ind=data(j).rate.*data(j).time>5e5 & isfinite(data(j).rate); data(j).auto1 = data(j).auto1(:,ind); data(j).auto2 = data(j).auto2(:,ind); data(j).tcspc1 = data(j).tcspc1(:,ind); data(j).tcspc2 = data(j).tcspc2(:,ind); data(j).time = data(j).time(ind); data(j).rate = data(j).rate(ind); end
% save data data hdata names
% for j=1:length(names) name{j}=names(j).name; end

for j=1:length(names) auto(:,j)=sum(data(j).auto1')' + sum(data(j).auto2')'; end
mm = length(sum(data(1).tcspc1,2));
for j=2:length(names)
    mm = min([mm length(sum(data(j).tcspc1,2))]);
end
for j=1:length(names)
    tcspc1(:,j) = sum(data(j).tcspc1(15:mm,:)')';
    tcspc2(:,j) = sum(data(j).tcspc2(15:mm,:)')';
end
tau = data(1).tau(15:mm);
resolution = hdata(1).Resolution;
ind = tau>117;
tau = tau(ind);
tcspc1 = tcspc1(ind,:);
tcspc2 = tcspc2(ind,:);
    
% DDM buffer
jbuf = 1;
jgfp = 7;
jp2x = 4;
jgly = 2;
jpro = [jgfp jp2x jgly];

pbuffer1 = simplex('expfun',1,0,[],[],[],tau,tcspc1(:,jbuf));
pbuffer2 = simplex('expfun',1,0,[],[],[],tau,tcspc2(:,jbuf));

[err, cbuffer1, zbuffer1] = expfun(pbuffer1,tau,tcspc1(:,jbuf));
[err, cbuffer2, zbuffer2] = expfun(pbuffer2,tau,tcspc2(:,jbuf));

clear p1 p2 c1 c2 tmp1 tmp2 z1 z2
for j=1:length(jpro)
    p1(:,j) = simplex('expfun',[1 pbuffer1],[0 pbuffer1],[inf pbuffer1],[],[],tau,tcspc1(:,jpro(j)));
    p2(:,j) = simplex('expfun',[1 pbuffer2],[0 pbuffer2],[inf pbuffer2],[],[],tau,tcspc2(:,jpro(j)));
    [err, c1(:,j), tmp1] = expfun(p1(:,j),tau,tcspc1(:,jpro(j)));
    [err, c2(:,j), tmp2] = expfun(p2(:,j),tau,tcspc2(:,jpro(j)));
    z1(:,j) = tmp1(:,1);
    z2(:,j) = tmp2(:,1);
end
% semilogy(tau,z1,tau,z2,tau,zbuffer1,tau,zbuffer2)

for j=1:length(jpro);
    para = [tau'/resolution tau'/resolution z1(:,j) z2(:,j) zbuffer1 zbuffer2 ones(size(zbuffer1)) ones(size(zbuffer2))];
    protein(j) = tttrpro(names(jpro(j)).name,'tcspcfcsfilter',para);
    save ddm protein
end

% PBS buffer
jbuf = 6;
jgfp = 8;
jp2x = 5;
jgly = 3;
jpro = [jgfp jp2x jgly];

pbuffer1 = simplex('expfun',1,0,[],[],[],tau,tcspc1(:,jbuf));
pbuffer2 = simplex('expfun',1,0,[],[],[],tau,tcspc2(:,jbuf));

[err, cbuffer1, zbuffer1] = expfun(pbuffer1,tau,tcspc1(:,jbuf));
[err, cbuffer2, zbuffer2] = expfun(pbuffer2,tau,tcspc2(:,jbuf));

clear p1 p2 c1 c2 tmp1 tmp2 z1 z2
for j=1:length(jpro)
    p1(:,j) = simplex('expfun',[1 pbuffer1],[0 pbuffer1],[inf pbuffer1],[],[],tau,tcspc1(:,jpro(j)));
    p2(:,j) = simplex('expfun',[1 pbuffer2],[0 pbuffer2],[inf pbuffer2],[],[],tau,tcspc2(:,jpro(j)));
    [err, c1(:,j), tmp1] = expfun(p1(:,j),tau,tcspc1(:,jpro(j)));
    [err, c2(:,j), tmp2] = expfun(p2(:,j),tau,tcspc2(:,jpro(j)));
    z1(:,j) = tmp1(:,1);
    z2(:,j) = tmp2(:,1);
end
% semilogy(tau,z1,tau,z2,tau,zbuffer1,tau,zbuffer2)

for j=1:length(jpro);
    para = [tau'/resolution tau'/resolution z1(:,j) z2(:,j) zbuffer1 zbuffer2 ones(size(zbuffer1)) ones(size(zbuffer2))];
    protein(j) = tttrpro(names(jpro(j)).name,'tcspcfcsfilter',para);
    save pbs protein
end



% AUSWERTUNG

load pbs
jbuf = 6;
jgfp = 8;
jp2x = 5;
jgly = 3;
jpro = [jgfp jp2x jgly];
t = protein(1).autotime;
for j=1:length(jpro)
    ind = protein(j).rate.*protein(j).time<2e5; 
    protein(j).auto(:,:,:,ind) = []; 
    protein(j).auto1(:,:,:,ind) = []; 
    protein(j).auto2(:,:,:,ind) = []; 
    auto0(:,:,:,j) = sum(protein(j).auto,4);
    auto1(:,:,:,j) = sum(protein(j).auto1,4);
    auto2(:,:,:,j) = sum(protein(j).auto2,4);    
end
% semilogx(t,auto0(:,1,1,1),t,auto0(:,2,2,1),t,auto0(:,3,3,1))


load ddm
jbuf = 1;
jgfp = 7;
jp2x = 4;
jgly = 2;
jpro = [jgfp jp2x jgly];
t = protein(1).autotime;
for j=1:length(jpro)
    ind = protein(j).rate.*protein(j).time<2e5; 
    protein(j).auto(:,:,:,ind) = []; 
    protein(j).auto1(:,:,:,ind) = []; 
    protein(j).auto2(:,:,:,ind) = []; 
    auto0(:,:,:,j) = sum(protein(j).auto,4);
    auto1(:,:,:,j) = sum(protein(j).auto1,4);
    auto2(:,:,:,j) = sum(protein(j).auto2,4);    
end
% semilogx(t,auto0(:,1,1,1),t,auto0(:,2,2,1),t,auto0(:,3,3,1))



clear p1 p2 c1 c2 
pguess = [2.5e-3 1e-4 2.5e5];
for j=1:length(jpro)
    p(:,j) = simplex('rigler',pguess,[0 0 0],[],[],[],t,auto(:,jpro(j)));
    [err c(:,j)] = rigler(p(:,j),t,auto(:,jpro(j)));
    p0(:,j) = simplex('rigler',p(:,j),[0 0 0],[],[],[],t,auto0(:,1,1,j));
    [err c0(:,j)] = rigler(p0(:,j),t,auto0(:,1,1,j));
    p1(:,j) = simplex('rigler',p(:,j),[0 0 0],[],[],[],t,auto1(:,1,1,j));
    [err c1(:,j)] = rigler(p1(:,j),t,auto1(:,1,1,j));
    p2(:,j) = simplex('rigler',p(:,j),[0 0 0],[],[],[],t,auto2(:,1,1,j));
    [err c2(:,j)] = rigler(p2(:,j),t,auto2(:,1,1,j));    
end

for k=1:length(jpro) 
    conc(k)=c(1,k)/c(2,k); 
    conc0(k)=c0(1,k)/c0(2,k); 
    conc1(k)=c1(1,k)/c1(2,k); 
    conc2(k)=c2(1,k)/c2(2,k); 
end;
for k=1:length(jpro) 
    bright(k)=sqrt(c(1,k))/conc(k);     
    bright0(k)=sqrt(c0(1,k))/conc0(k); 
    bright1(k)=sqrt(c1(1,k))/conc1(k); 
    bright2(k)=sqrt(c2(1,k))/conc2(k);     
end;


