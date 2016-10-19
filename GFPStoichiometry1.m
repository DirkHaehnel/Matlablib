% [buffer hbuffer]=tttrpro('E:\041111\Puffer_01.t3r','xfcsfilter');
% [gfp hgfp]=tttrpro('E:\041111\egfp_01.t3r','xfcsfilter');
% [p2x1 hp2x1]=tttrpro('E:\041111\P2X1_01.t3r','xfcsfilter');
% [glyr hglyr]=tttrpro('E:\041111\GlyR_01.t3r','xfcsfilter');
% save gfpStoichiometry buffer hbuffer gfp hgfp p2x1 hp2x1 glyr hglyr

auto = [sum(buffer.auto1(:,1:17),2)+sum(buffer.auto2(:,1:17),2), sum(gfp.auto1(:,1:end-2),2)+sum(gfp.auto2(:,1:end-2),2), ...
        sum(glyr.auto1(:,1:40),2)+sum(glyr.auto2(:,1:40),2), sum(p2x1.auto1(:,1:18),2)+sum(p2x1.auto2(:,1:18),2)];
mm = min([length(sum(buffer.tcspc1(:,1:17),2)), length(sum(gfp.tcspc1(:,1:end-2),2)), length(sum(glyr.tcspc1(:,1:40),2)), length(sum(p2x1.tcspc1(:,1:18),2))]);
tcspc1 = [sum(buffer.tcspc1(15:mm,1:17),2), sum(gfp.tcspc1(15:mm,1:end-2),2), sum(glyr.tcspc1(15:mm,1:40),2), sum(p2x1.tcspc1(15:mm,1:18),2)];
tcspc2 = [sum(buffer.tcspc2(15:mm,1:17),2), sum(gfp.tcspc2(15:mm,1:end-2),2), sum(glyr.tcspc2(15:mm,1:40),2), sum(p2x1.tcspc2(15:mm,1:18),2)];
tau = buffer.tau(15:mm);

ind = tau>118;
pbuffer1 = simplex('expfun',1./[1 6],[0 0],[],[],[],tau(ind),tcspc1(ind,1));
pbuffer2 = simplex('expfun',1./[1 6],[0 0],[],[],[],tau(ind),tcspc2(ind,1));
for k=1:3
    pbuffer1 = simplex('expfun',pbuffer1,0*pbuffer1,[],[],[],tau(ind),tcspc1(ind,1));
    pbuffer2 = simplex('expfun',pbuffer2,0*pbuffer1,[],[],[],tau(ind),tcspc2(ind,1));
end
[err, cbuffer1, zbuffer1] = expfun(pbuffer1,tau(ind),tcspc1(ind,1));
[err, cbuffer2, zbuffer2] = expfun(pbuffer2,tau(ind),tcspc2(ind,1));
zbuffer1 = zbuffer1*cbuffer1(2:end);
zbuffer2 = zbuffer2*cbuffer2(2:end);

clear p1 p2 c1 c2 tmp1 tmp2 z1 z2
for j=2:4
    p1(:,j-1)=simplex('expfun',[0.5 0.1],[0 0],[],[],[],tau(ind),tcspc1(ind,j));
    p2(:,j-1)=simplex('expfun',[0.5,0.1],[0 0],[],[],[],tau(ind),tcspc2(ind,j));
    [err, c1(:,j-1), tmp1] = expfun(p1(:,j-1),tau(ind),tcspc1(ind,j));
    [err, c2(:,j-1), tmp2] = expfun(p2(:,j-1),tau(ind),tcspc2(ind,j));
    z1(:,2*j-3:2*j-2) = tmp1;
    z2(:,2*j-3:2*j-2) = tmp2;
end

j = 1;
para = [tau(ind)'/hbuffer.Resolution tau(ind)'/hbuffer.Resolution z1(:,2*j-1) z2(:,2*j-1) z1(:,2*j) z2(:,2*j) ones(size(zbuffer1)) ones(size(zbuffer2))];
gfptcspc=tttrpro('f:\041111\egfp_01.t3r','tcspcfcsfilter',para);
j = 2;
para = [tau(ind)'/hbuffer.Resolution tau(ind)'/hbuffer.Resolution z1(:,2*j-1) z2(:,2*j-1) z1(:,2*j) z2(:,2*j) ones(size(zbuffer1)) ones(size(zbuffer2))];
glyrtcspc=tttrpro('f:\041111\GlyR_01.t3r','tcspcfcsfilter',para);
j = 3;
para = [tau(ind)'/hbuffer.Resolution tau(ind)'/hbuffer.Resolution z1(:,2*j-1) z2(:,2*j-1) z1(:,2*j) z2(:,2*j) ones(size(zbuffer1)) ones(size(zbuffer2))];
p2x1tcspc=tttrpro('f:\041111\P2X1_01.t3r','tcspcfcsfilter',para);

gfp11 = [sum(gfptcspc.auto(:,1,1,1:end-2),4) sum(gfptcspc.auto1(:,1,1,1:end-2),4) sum(gfptcspc.auto2(:,1,1,1:end-2),4)];
glyr11 = [sum(glyrtcspc.auto(:,1,1,1:45),4) sum(glyrtcspc.auto1(:,1,1,1:45),4) sum(glyrtcspc.auto2(:,1,1,1:45),4)];
p2x111 = [sum(p2x1tcspc.auto(:,1,1,1:18),4) sum(p2x1tcspc.auto1(:,1,1,1:18),4) sum(p2x1tcspc.auto2(:,1,1,1:18),4)];

t = gfptcspc.autotime;
auto  = [sum(gfp.auto1(:,1:end-2),2)+sum(gfp.auto2(:,1:end-2),2) sum(glyr.auto1(:,1:40),2)+sum(glyr.auto2(:,1:40),2), ...
        sum(p2x1.auto1(:,1:18),2)+sum(p2x1.auto2(:,1:18),2)];

p(:,1,1) = simplex('rigler',[1e-2 1e-3],[0 0],[],[],[],t(6:end),gfp11(6:end,1));
p(:,2,1) = simplex('rigler',[1e-2 1e-3],[0 0],[],[],[],t(6:end),gfp11(6:end,2));
p(:,3,1) = simplex('rigler',[1e-2 1e-3],[0 0],[],[],[],t(6:end),gfp11(6:end,3));
[err c(:,1,1)] = rigler(p(:,1,1),t(6:end),gfp11(6:end,1));
[err c(:,2,1)] = rigler(p(:,2,1),t(6:end),gfp11(6:end,2));
[err c(:,3,1)] = rigler(p(:,3,1),t(6:end),gfp11(6:end,3));

p(:,1,2) = simplex('rigler',[1e-2 1e-3],[0 0],[],[],[],t(6:end),p2x111(6:end,1));
p(:,2,2) = simplex('rigler',[1e-2 1e-3],[0 0],[],[],[],t(6:end),p2x111(6:end,2));
p(:,3,2) = simplex('rigler',[1e-2 1e-3],[0 0],[],[],[],t(6:end),p2x111(6:end,3));
[err c(:,1,2)] = rigler(p(:,1,2),t(6:end),p2x111(6:end,1));
[err c(:,2,2)] = rigler(p(:,2,2),t(6:end),p2x111(6:end,2));
[err c(:,3,2)] = rigler(p(:,3,2),t(6:end),p2x111(6:end,3));

p(:,1,3) = simplex('rigler',[1e-2 1e-3],[0 0],[],[],[],t(6:end),glyr11(6:end,1));
p(:,2,3) = simplex('rigler',[1e-2 1e-3],[0 0],[],[],[],t(6:end),glyr11(6:end,2));
p(:,3,3) = simplex('rigler',[1e-2 1e-3],[0 0],[],[],[],t(6:end),glyr11(6:end,3));
[err c(:,1,3)] = rigler(p(:,1,3),t(6:end),glyr11(6:end,1));
[err c(:,2,3)] = rigler(p(:,2,3),t(6:end),glyr11(6:end,2));
[err c(:,3,3)] = rigler(p(:,3,3),t(6:end),glyr11(6:end,3));

for j=1:3 for k=1:3 conc(j,k)=c(1,j,k)/c(2,j,k); end; end
for j=1:3 for k=1:3 bright(j,k)=sqrt(c(1,j,k))/conc(j,k); end; end



% [err c(:,1,1)] = rigler(p0(:,1),t(6:end),gfp11(6:end,1));
% [err c(:,2,1)] = rigler(p0(:,1),t(6:end),gfp11(6:end,2));
% [err c(:,3,1)] = rigler(p0(:,1),t(6:end),gfp11(6:end,3));
% [err c(:,1,2)] = rigler(p0(:,2),t(6:end),p2x111(6:end,1));
% [err c(:,2,2)] = rigler(p0(:,2),t(6:end),p2x111(6:end,2));
% [err c(:,3,2)] = rigler(p0(:,2),t(6:end),p2x111(6:end,3));
% [err c(:,1,3)] = rigler(p0(:,3),t(6:end),glyr11(6:end,1));
% [err c(:,2,3)] = rigler(p0(:,3),t(6:end),glyr11(6:end,2));
% [err c(:,3,3)] = rigler(p0(:,3),t(6:end),glyr11(6:end,3));
% 
% for j=1:3 for k=1:3 conc(j,k)=c(1,j,k)/c(2,j,k); end; end
% for j=1:3 for k=1:3 bright(j,k)=sqrt(c(1,j,k))/conc(j,k); end; end


% new try


% j = 1;
% para = [tau(ind)'/hbuffer.Resolution tau(ind)'/hbuffer.Resolution ind(ind)' ind(ind)'];
% gfpfilt=tttrpro('E:\041111\egfp_01.t3r','tcspcfcsfilter',para);
% save gfpStoichiometry gfpfilt -append
% j = 2;
% glyrfilt=tttrpro('E:\041111\GlyR_01.t3r','tcspcfcsfilter',para);
% save gfpStoichiometry glyrfilt -append
% j = 3;
% p2x1filt=tttrpro('E:\041111\P2X1_01.t3r','tcspcfcsfilter',para);
% save gfpStoichiometry p2x1filt -append
% 
% t = gfpfilt.autotime;
% auto  = [sum(gfpfilt.auto(:,1,1,1:end-2),4) sum(glyrfilt.auto(:,1,1,1:40),4) sum(p2x1filt.auto(:,1,1,1:18),4)];
% 
% auto1  = [sum(gfpfilt.auto1(:,1,1,1:end-2),4) sum(glyrfilt.auto1(:,1,1,1:40),4) sum(p2x1filt.auto1(:,1,1,1:18),4)];
% pguess = [2.5e-3 1e-4 2.5e5];
% p0(:,1) = simplex('rigler',pguess,[0 0 0],[],[],[],t(6:end),auto(6:end,1));
% [err c0(:,1)] = rigler(p0(:,1),t(6:end),auto(6:end,1));
% p0(:,2) = simplex('rigler',pguess,[0 0 0],[],[],[],t(6:end),auto(6:end,2));
% [err c0(:,2)] = rigler(p0(:,2),t(6:end),auto(6:end,2));
% p0(:,3) = simplex('rigler',pguess,[0 0 0],[],[],[],t(6:end),auto(6:end,3));
% [err c0(:,3)] = rigler(p0(:,3),t(6:end),auto(6:end,3));
% pbuffer=simplex('rigler',[1e-2 1e-3],[0 0],[],[],[],buffer.autotime(6:end),sum(buffer.auto1(6:end,1:17),2)+sum(buffer.auto2(6:end,1:17),2));
% [err cbuffer]=rigler(pbuffer,buffer.autotime(6:end),sum(buffer.auto1(6:end,1:17),2)+sum(buffer.auto2(6:end,1:17),2));
% for k=1:3 tmp(k)=(c1(3,k)/sum(c1(:,k)))^2*c0(1,k)/cbuffer(1)*cbuffer(2); end;
% for k=1:3 c0(1,k)=(c1(2,k)/sum(c1(:,k)))^2*c0(1,k); end;
% for k=1:3 c0(2,k)=c0(2,k)-tmp(k); end;
% for k=1:3 conc0(k)=c0(1,k)/c0(2,k); end;
% for k=1:3 bright0(k)=sqrt(c0(1,k))/conc0(k); end;
% 
% auto2  = [sum(gfpfilt.auto2(:,1,1,1:end-2),4) sum(glyrfilt.auto2(:,1,1,1:40),4) sum(p2x1filt.auto2(:,1,1,1:18),4)];
% pguess = [2.5e-3 1e-4 2.5e5];
% p0(:,1) = simplex('rigler',pguess,[0 0 0],[],[],[],t(6:end),auto(6:end,1));
% [err c0(:,1)] = rigler(p0(:,1),t(6:end),auto(6:end,1));
% p0(:,2) = simplex('rigler',pguess,[0 0 0],[],[],[],t(6:end),auto(6:end,2));
% [err c0(:,2)] = rigler(p0(:,2),t(6:end),auto(6:end,2));
% p0(:,3) = simplex('rigler',pguess,[0 0 0],[],[],[],t(6:end),auto(6:end,3));
% [err c0(:,3)] = rigler(p0(:,3),t(6:end),auto(6:end,3));
% pbuffer=simplex('rigler',[1e-2 1e-3],[0 0],[],[],[],buffer.autotime(6:end),sum(buffer.auto1(6:end,1:17),2)+sum(buffer.auto2(6:end,1:17),2));
% [err cbuffer]=rigler(pbuffer,buffer.autotime(6:end),sum(buffer.auto1(6:end,1:17),2)+sum(buffer.auto2(6:end,1:17),2));
% for k=1:3 tmp(k)=(c2(3,k)/sum(c2(:,k)))^2*c0(1,k)/cbuffer(1)*cbuffer(2); end;
% for k=1:3 c0(1,k)=(c2(2,k)/sum(c2(:,k)))^2*c0(1,k); end;
% for k=1:3 c0(2,k)=c0(2,k)-tmp(k); end;
% for k=1:3 conc0(k)=c0(1,k)/c0(2,k); end;
% for k=1:3 bright0(k)=sqrt(c0(1,k))/conc0(k); end;


