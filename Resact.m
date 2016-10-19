% program for analyzing Resact dynamics

load D:\Daten\Gregor\Resact\Resact.mat
NatConst

Temp = 271.15 + TempW9;
t = resW9(1).autotime(2:end);
len = length(Temp);

for j=1:len
    auto(:,j) = mean(resW9(j).auto(2:end,:),2);
end

ind = t<5e-6 | t>1;
for j=1:len
    p(:,j) = Simplex('rigler',[1e-3 1e-4 1e6],[],[],[],[],t(ind),auto(ind,j));
    for k=1:3
        p(:,j) = Simplex('rigler',p(:,j),[],[],[],[],t(ind),auto(ind,j));
    end
    p(:,j) = Simplex('rigler',p(:,j),[],[],[],[],t(ind),auto(ind,j),[],1);
    [err, c(:,j)] = rigler(p(:,j),t(ind),auto(ind,j));
    z(:,j) = rigler(p(:,j),c(:,j),t);
end

semilogx(t,auto,'o',t,z)

kp = p(3,:).*c(3,:)./sum(c(2:3,:));
km = p(3,:).*c(2,:)./sum(c(2:3,:));

% ind = [1 2 3 5 6];
ind = 1:6;

ap = polyfit(1./Temp(ind),log(kp(ind)),1);
am = polyfit(1./Temp(ind),log(km(ind)),1);

plot(1./Temp(ind),kp(ind),'-o',1./Temp(ind),km(ind),'-o',1./Temp,exp(polyval(ap,1./Temp)),1./Temp,exp(polyval(am,1./Temp)))
legend({['\itk_p\rm = ' mnum2str(-ap(1)*MolarGasConstant/1e3,2,2) ' kJ/mol'],['\itk_m\rm = ' mnum2str(-am(1)*MolarGasConstant/1e3,2,2)  ' kJ/mol']})
xlabel('1/\itT\rm [1/K]')
ylabel('rate constant')

return

pp = Simplex('ResactFun',[reshape(p(1:2,:),1,2*len) exp(ap(2)) exp(am(2)) -ap(1) -am(1)],[],[],[],[],t,auto,Temp,1)

