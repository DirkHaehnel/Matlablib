clear all
load 'Y:\GCAT\lt und wichtung.mat'

tauval = [2 4];

tau = (1:length(tcspc))';
plot(tau,tcspc)
[a,bla] = ginput(2);

ind = tau>a(1) & tau<a(2);
tau0 = tau(ind);
tau1 = tau0(1)-1;

master = exp(-head.Resolution*(tau0-tau1)*(1./tauval));
master = master./(ones(size(master,1),1)*sum(master));

p1 = zeros(size(lt));
p2 = p1;
for j=1:size(lt,1) 
    for k=1:size(lt,2) 
        tmp = lt{j,k};
        if ~isempty(tmp)
            tmp(tmp<tau0(1) | tmp>tau0(end)) = [];
        end
        if ~isempty(tmp)
            p1(j,k) = sum(log(master(tmp-tau1,1)))/length(tmp); 
            p2(j,k) = sum(log(master(tmp-tau1,2)))/length(tmp);             
        end
    end; 
end

tmp = sum(p1-p2);
tmp = (tmp-min(tmp))/(max(tmp)-min(tmp));
plot(1:size(lt,2),tmp)
xlabel('position along DNA');
ylabel('short lifetime / long lifetime')


