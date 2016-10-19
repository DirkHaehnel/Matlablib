clear all
fname = 'D:\Joerg\Doc\Defocused\Hiroshi\t27_5_003M001_res.mat'; % file name of data file
eval(['load ' fname ' res'])
kv = 1:50; % max number of frame number difference to be analyzed
clear pp
maxord = 4;

%%%% following line eliminates non-valid results %%%
res(res(:,4)==0 & res(:,5)==0,:) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1
    tst = ones(size(res,1),1);
    tst(abs(diff(res([1 1:end],5)))>180)=-1;
    tst = cumprod(tst);
    theta = (90*(1-tst)+res(:,4).*tst)/180*pi;
    theta = unwrap(2*theta)/2;
else
    theta = res(:,4)/180*pi;
end
phi = unwrap(res(:,5)/180*2*pi)/2;
t = 1:size(res,1);
t = t(~(res([1:end-1],1)+1==res([2:end],1)));
t = [0 t size(res,1)];
for j=1:length(t)-1
    pp(j).theta = theta(t(j)+1:t(j+1));
    pp(j).phi = phi(t(j)+1:t(j+1));
end

j = 1;
tmp = [sin(pp(j).theta).*cos(pp(j).phi), sin(pp(j).theta).*sin(pp(j).phi), cos(pp(j).theta)];
for k=1:max(kv)
    if k<size(tmp,1)
        cnt(k) = size(tmp,1)-k;
        for s=1:maxord
            mm(k,s) = sum(sum((tmp(1:end-k,:).*tmp(k+1:end,:)),2).^s);
        end
    else
        for s=1:maxord
            mm(k,s) = 0;
        end
        cnt(k) = 0;
    end
end
for j=2:length(pp)
    tmp = [sin(pp(j).theta).*cos(pp(j).phi), sin(pp(j).theta).*sin(pp(j).phi), cos(pp(j).theta)];
    for k=1:max(kv)
        if k<size(tmp,1)
            cnt(k) = cnt(k) + size(tmp,1)-k;
            for s=1:maxord
                mm(k,s) = mm(k,s) + sum(sum((tmp(1:end-k,:).*tmp(k+1:end,:)),2).^s);
            end
        end
    end
end
mm = mm./(cnt'*ones(1,maxord));

% para = [1e-1 1e-2 1e-3 1e-4];
para = [1e-1 1e-2];
maxk = 25;
close all
for j=1:5
    para = simplex('RotoDiffFit',para,0*para,[],[],[],[kv(1:maxk)],[mm(1:maxk,:)]),
end
h=get(gca,'children');
for j=1:4
    set(h(j),'color',get(h(j+4),'color'));
end
xlabel('frame number'); ylabel('\langlecos^{\itn\rm}\Theta\rangle');
legend({'\langlecos\Theta\rangle','\langlecos^2\Theta\rangle','\langlecos^3\Theta\rangle','\langlecos^4\Theta\rangle'},3)

% plot sphere plot
figure
sphere(50)
alpha(0.1)
shading flat
hold on
j = 1;
tmp = [sin(pp(j).theta).*cos(pp(j).phi), sin(pp(j).theta).*sin(pp(j).phi), cos(pp(j).theta)];
plot3(tmp(:,1),tmp(:,2),tmp(:,3));
for j=2:length(pp)
    tmp = [sin(pp(j).theta).*cos(pp(j).phi), sin(pp(j).theta).*sin(pp(j).phi), cos(pp(j).theta)];
    plot3(tmp(:,1),tmp(:,2),tmp(:,3));
end
hold off
axis equal
lighting gouraud
axis off
h = line([0 0 0; 1.3 0 0],[0 0 0; 0 1.3 0],[0 0 0; 0 0 1.3]);
for j=1:length(h) set(h,'color','b','linewidth',1.5); end
view([135 30])




return % make a movie

for jj=1:36
    phi0 = (jj-1)*pi/18;
    sphere(50)
    alpha(0.1)
    shading flat
    hold on
    j = 1;
    tmp = [sin(pp(j).theta).*cos(pp(j).phi-phi0), sin(pp(j).theta).*sin(pp(j).phi-phi0), cos(pp(j).theta)];
    plot3(tmp(:,1),tmp(:,2),tmp(:,3),'linewidth',1);
    for j=2:length(pp)
        tmp = [sin(pp(j).theta).*cos(pp(j).phi-phi0), sin(pp(j).theta).*sin(pp(j).phi-phi0), cos(pp(j).theta)];
        plot3(tmp(:,1),tmp(:,2),tmp(:,3),'linewidth',1);
    end
    hold off
    axis equal
    lighting gouraud
    axis off
    h = line([0 0 0; 1.3*cos(phi0) 1.3*sin(phi0) 0],[0 0 0; -1.3*sin(phi0) 1.3*cos(phi0) 0],[0 0 0; 0 0 1.3]);
    for j=1:length(h) set(h,'color','b','linewidth',1.5); end
    view([135 30])
    eval(['print -dpng tmp' mint2str(jj,2)])
end



