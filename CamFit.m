function [err, c, z] = CamFit(p, conc, y, conc0, bld)

conc = conc(:);
y = y(:);

z = zeros(length(conc),length(p));
for j=1:length(p)
    z(:,j) = 1./(p(j)./conc+1).^2;
end
c = [ones(size(conc)) z]\y;
z = [ones(size(conc)) z]*c;

if nargin>4 && bld==1
    semilogx(conc, y, 'o', conc, z); drawnow;
end

err = abs(sum((y-z).^2./abs(z)));

if nargin>3 && ~isempty(conc0)
    conc0 = conc0(:);
    z = zeros(length(conc0),length(p));
    for j=1:length(p)
        z(:,j) = 1./(p(j)./conc0+1).^2;
    end
    z = [ones(size(conc0)) z]*c;
    semilogx(conc, y, 'o', conc0, z); drawnow;
end

    
return

p=Simplex('CamFit',[1e-6 1e-4],[0 0],[],[],[], c, rh);
errorbar(c,rh/10,drh/10,'o'); hold on; plot(10.^(-8:0.01:-2),CamFit(p,c,rh/10,10.^(-8:0.01:-2))); hold off
set(gca,'xscale','log')
axis([1e-8 1e-2 2.15 2.45])
xlabel('Ca^{2+} concentration (M)')
ylabel('hydrodynamic radius (nm)')
