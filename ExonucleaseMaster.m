% Exonucelase Master

t = 0:2e-3:1;
k0 = 108; sig = 66;
kv = 10:300;
weight = exp(-(kv-k0).^2/2/sig^2);
seq = [1, 2, 2, 1, 4, 1, 2, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 2, 1, 2, 2, 2, 2, 1, 1, 4, 3];
nv = 1:length(seq);
nv = nv(seq==max(seq));

z = zeros(length(seq),length(t)); 
for j=1:length(kv)
    z = z + weight(j)*ExonucleaseFun(kv(j), t, ones(size(seq)), [], nv);
end


close all; p0 = simplex('ExonucleaseFun',100,[0],[300],[],[],t(1:200),ones(size(seq)),sum(z(nv,1:200)),nv);
[err,z0,z0] = ExonucleaseFun(p0, t, ones(size(seq)), sum(z(nv,:)), nv);

tmp = seq; tmp(tmp==2)=1; tmp(tmp==3)=1; tmp(tmp==4)=2;
close all; p1 = simplex('ExonucleaseFun',[300, 300; 10, 1],[],[],[],[],t,tmp,sum(z(nv,:)));
close all; for j=1:5 p1 = simplex('ExonucleaseFun',p1,0*p1,[],[],[],t,tmp,sum(z(nv,:))); end
[err,z1,z1] = ExonucleaseFun(p1, t, tmp, sum(z(nv,:)));

tmp = seq; tmp(tmp==3)=1; tmp(tmp==4)=3;
close all; p1a = simplex('ExonucleaseFun',[300, 300, 10],[],[],[],[],t,tmp,sum(z(nv,:)));
close all; for j=1:5 p1a = simplex('ExonucleaseFun',p1a,0*p1a,[],[],[],t,tmp,sum(z(nv,:))); end
[err,z1a,z1a] = ExonucleaseFun(p1a, t, tmp, sum(z(nv,:)));

close all; p2 = simplex('ExonucleaseFun',[300, 300 300 300; 300 300 300 300; 300, 300, 300, 10; 10, 10, 10, 1],[],[],[],[],t,seq,sum(z(nv,:)));
close all; for j=1:5 p2 = simplex('ExonucleaseFun',p2,0*p2,[],[],[],t,seq,sum(z(nv,:))); end
p2 = reshape(p2,4,4);
[err,z2,z2] = ExonucleaseFun(p2, t, seq, sum(z(nv,:)));

save Exonuclease p0 p1 p2 z seq nv

% Fig2
plot(t,sum(z(nv,:)),'o',t,z0,t,z1,t,z2)

return

% Fig3
tmp1 = seq; tmp1(7:36) = 1;
tmp2 = seq; tmp2(7:36) = 2;
[z2a, M, weight] = ExonucleaseFun(p2, t, seq);
z2a = weight*z2a(nv,:);
z2b = ExonucleaseFun(p2, t, tmp1);
z2b = weight*z2b(nv,:);
z2c = ExonucleaseFun(p2, t, tmp2);
z2c = weight*z2c(nv,:);

plot(t,z2a,t,z2b,t,z2c)
