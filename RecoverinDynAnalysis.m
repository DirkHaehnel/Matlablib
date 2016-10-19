cd E:\Daten\Gregor\Recoverin\
names = dir('*.t3r');

for j=1:length(names)
    [recodyn(j), hrecodyn(j)] = Cy5dynamics(names(j).name);
end

save recodyndata recodyn hrecodyn

return


t = cydyn(1).autotime;
ind = t>=1e-6;
t = t(ind);
for j=1:length(cydyn)
    x(:,:,:,j) = cydyn(j).auto(ind,1:2,1:2);
end

for j=1:size(x,4)
    err = 0;
    tmp = [1e-2 1e-3 1e-6 100e-6 10e-6]';
    for casc=1:10
        p0 = tmp(:,err==min(err));
        for sub=1:10
            tmp(:,sub) = simplex('FLCS',p0.*exp(log(2)*(2*rand(size(p0))-1)),zeros(size(p0)),[],[],[],t,x(:,:,:,j));
            err(sub) = FLCS(tmp(:,sub),t,x(:,:,:,j)); 
        end
    end
    p(:,j) = tmp(:,err==min(err));
    [err, c(:,:,j)] = FLCS(p(:,j),t,x(:,:,:,j));
end

% save Cy5DynAnalsyis p c name

plot(pow([1 2 5]),1./p(5,[1 2 5]),'o',pow([1 2 5]),polyval(polyfit(pow([1 2 5]),1./p(5,[1 2 5]),1),pow([1 2 5])))
xlabel('power [\muW]')
ylabel('transition rate \itk_{sn}\rm + \itk_{ns}\rm [1/s]')

for j=1:5 r(j) = sqrt(c(5,2,j)/c(5,3,j)/c(5,4,j)*c(5,1,j)); end
plot(pow([1 2 5]),r([1 2 5]),'o',0:400,ones(size(0:400))*mean(r([1 2 5])))
axis([0 400 0 4])
xlabel('power [\muW]')
ylabel('transition rate ratio \itk_{sn}\rm/\itk_{ns}\rm')

ksn = r([1 2 5])./(1+r([1 2 5]))./p(5,[1 2 5]);
kns = 1./(1+r([1 2 5]))./p(5,[1 2 5]);
ksn400 = polyval(polyfit(pow([1 2 5]),ksn,1),0:400);
kns400 = polyval(polyfit(pow([1 2 5]),kns,1),0:400);
ksn400a = polyval(polyfit([-pow([1 2 5]) pow([1 2 5])],[-ksn ksn],1),0:400);
kns400a = polyval(polyfit([-pow([1 2 5]) pow([1 2 5])],[-kns kns],1),0:400);

plot(0:400,ksn400,'r',0:400,kns400,'b',pow([1 2 5]),ksn,'or',pow([1 2 5]),kns,'ob')
xlabel('power [\muW]')
ylabel('transition rate [1/s]')
legend({'\itk_{sn}','\itk_{ns}'},2)

plot(0:400,ksn400,'r',0:400,kns400,'b',0:400,ksn400a,'r:',0:400,kns400a,'b:',pow([1 2 5]),ksn,'or',pow([1 2 5]),kns,'ob')
legend({'\itk_{sn}','\itk_{ns}'},2)
xlabel('power [\muW]')
ylabel('transition rate [1/s]')

