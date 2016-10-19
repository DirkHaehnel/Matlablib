% Program for Daniel Koller: Lifetime calculation in layered structures

clear all
close all

n1 = 1.515;
nv = [[1.5707 1.5713 1.578]; [1.5664 1.567 1.5735]];
n2 = 1;

dv = [1 10];
lamemv = [560 610]/1e3;

zv = (0.5:3e2)/3e2; col = ones(size(zv));

for jd = 2:length(dv)
    d = dv(jd);
    for jlam = 2:length(lamemv)
        lamem = lamemv(jlam);
        for jn = 2:2%size(nv,2)
            n = nv(jlam,jn);
            
            [lvd,lvu,lpd,lpu,qvd,qvu,qpd,qpu,qv,qp] = LifetimeL((zv+d-1)/lamem*2*pi,n1,n,n2,[],d/lamem*2*pi,[]);

            [lvd0,lvu0,lpd0,lpu0,qvd0,qvu0,qpd0,qpu0,qv0,qp0] = LifetimeL(0,[n1,n],n2,n2,d/lamem*2*pi,0,[]);
            
            plot(zv+d-1,4/3*n2./(qpu+qpd+qp),zv+d-1,4/3*n2./(qvu+qvd+qv),...
                zv+d-1,4/3*n2./(qpu0+qpd0+qp0)*col,':r',...
                zv+d-1,4/3*n2./(qvu0+qvd0+qv0)*col,':b')
            legend({'in-plane dipole','vertical dipole'},2)
            ylabel('lifetime (\tau_{\itair\rm})')
            xlabel('vertical position (\mum)')

            return
            
            eval(['print -dpng -r300 Lifetime_' mint2str(d,2) '_' mint2str(1e3*lamem,3) '_' mint2str(jn,1)])
            
            plot(zv+d-1,qpu./(qpu+qpd+qp),zv+d-1,qpd./(qpu+qpd+qp),zv+d-1,qp./(qpu+qpd+qp),'g',...
                zv+d-1,qpu0./(qpu0+qpd0+qp0)*col,':r',zv+d-1,qpd0./(qpu0+qpd0+qp0)*col,':b',zv+d-1,qp0./(qpu0+qpd0+qp0)*col,':g')
            legend({'emission into air','emission into substrate','emmission into guided modes'},2)
            ylabel('rel. emission')
            xlabel('vertical position (\mum)')

            eval(['print -dpng -r300 EmissionInPlane_' mint2str(d,2) '_' mint2str(1e3*lamem,3) '_' mint2str(jn,1)])            

            plot(zv+d-1,qvu./(qvu+qvd+qv),zv+d-1,qvd./(qvu+qvd+qv),zv+d-1,qv./(qvu+qvd+qv),'g',...
                zv+d-1,qvu0./(qvu0+qvd0+qv0)*col,':r',zv+d-1,qvd0./(qvu0+qvd0+qv0)*col,':b',zv+d-1,qv0./(qvu0+qvd0+qv0)*col,':g')
            legend({'emission into air','emission into substrate','emmission into guided modes'},2)
            ylabel('rel. emission')
            xlabel('vertical position (\mum)')

            eval(['print -dpng -r300 EmissionVertical_' mint2str(d,2) '_' mint2str(1e3*lamem,3) '_' mint2str(jn,1)])            
            
            close all
        end
    end
end

