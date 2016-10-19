% PicoHarpDeadTime
clear all

tau = 6; % decay time
trep = 12.5; % repetition period
td = 98; % dead time
dt = 0.01; % channel width

t = -td:dt:trep-dt;

pow = [0.01 0.1 1]/trep;
int0 = exp(-mod(t,trep)./tau)'/tau;

for jp = 1:length(pow)
    int(:,jp) = pow(jp)*int0;
    tt = t(t>0);
    for j=1:2
        for k=1:length(tt)
            int(t==tt(k),jp) = exp(-sum(int(t>=tt(k)-td & t<tt(k),jp))*dt)*int0(t==tt(k))*pow(jp);
        end
        for k=1:length(t<=0)
            int(k,jp) = int(round(t/dt)==round(mod(t(k),trep)/dt),jp);    
        end
        plot(t(t>0),int0(t>0),t(t>0),int(t>0,jp)/pow(jp)); drawnow
    end
end