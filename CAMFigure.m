load c:\Joerg\Doc\Fcs\2Focus\CAM\final_results.mat
[c,ind]=sort(c);rh=rh(ind);drh=drh(ind);
semilogx(c_fit,rh_fit)
hold on
errorbar(c,rh,drh,'ob')
hold off
axis([5e-9 1e-2 21.5 24.5])
xlabel('free Ca^{2+} concentration (M)')
ylabel(['hydrodynamic radius ' char(197)])