% figures for Atto655 paper

load D:\Joerg\Doc\Fcs\2Focus\Atto655\AttoData\GDCl_results.mat

c = polyfit(1./viscosity,twofocus_d,1);

v = [1./viscosity 1./nmr_visc];

plot(v,polyval(c,v),'b',1./viscosity,onefocus_d,'sg',1./viscosity,twofocus_d,'or',1./nmr_visc,nmr_d,'vc',...
    [1;1]*(1./viscosity),[onefocus_d-onefocus_std;onefocus_d+onefocus_std],'g',...
    [1;1]*(1./viscosity),[twofocus_d-twofocus_std;twofocus_d+twofocus_std],'r',...
    [1;1]./nmr_visc,[nmr_d*0.97;nmr_d*1.03],'c')
axis([0.6 1.7 1e-6 7e-6])   
grid
xlabel('1/viscosity [1/mPa\cdots]')
ylabel('diffusion coefficient [cm^2/s]')
legend({'lin fit','1-focus FCS','2-focus FCS','NMR'},4)

ref = polyval(c,1/viscosity(1))

load D:\Joerg\Doc\Fcs\2Focus\Atto655\AttoData\Opt_Sat_results.mat

c = polyfit(pow(pow<=30),twofocus_d(pow<=30),0);

plot([0 72.5],polyval(c,[0 72.5]),'b',pow,onefocus_d,'sg',pow,twofocus_d,'or',...
    [1;1]*pow,[onefocus_d-onefocus_std;onefocus_d+onefocus_std],'g',...
    [1;1]*pow,[twofocus_d-twofocus_std;twofocus_d+twofocus_std],'r')
axis([0 72.5 3.4e-6 4.8e-6])   
grid
xlabel('excitation power [\muW]')
ylabel('diffusion coefficient [cm^2/s]')
legend({'average','1-focus FCS','2-focus FCS'},2)	
