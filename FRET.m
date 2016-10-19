function [foerster, ab, em, ind] = FRET(em,ab,qf,ex,wl,RefractiveIndex,kappa2)

% FRET(em,ab,qf,ex,wl,RefractiveIndex) calculates the Foerster radius in nm
%
% em - emission spectrum of donor
% ab - absorption spectrum of acceptor
% qf - quantum yield of donor
% ex - molar extinction of acceptor (in l/mol/cm) at wavelength wl (in nm)

AvogadroConstant = 6.0221367*10.^23;
if nargin<6 || isempty(RefractiveIndex)
    RefractiveIndex = 1.33;
end

if nargin<7
    kappa2 = 2/3;
end

delta = 0.1;
ind = max([min(ab(:,1)) min(em(:,1))]):delta:min([max(ab(:,1)) max(em(:,1))]);
ab(:,2) = ab(:,2)/interp1(ab(:,1),ab(:,2),wl,'cubic')*ex*1e3*1e14*log(10)/AvogadroConstant; 
% 1e14 for converting cm^2 into nm^2
em = interp1(em(:,1),em(:,2),ind,'cubic');  

ab = interp1(ab(:,1),ab(:,2),ind,'cubic');

foerster = sum(ab.*em.*ind.^4)/sum(em);
foerster = kappa2*9/8/pi/((2*pi*RefractiveIndex)^4)*qf*foerster;
foerster = foerster^(1/6);

plot(ind,ab/max(ab),ind,em/max(em)); 
xlabel('wavelength [nm]'); 
ylabel('spectrum [a.u.]');
ax = axis;
text(ax(1)+0.6*(ax(2)-ax(1)),ax(3)+0.9*(ax(4)-ax(3)),['R_0 = ' num2str(foerster) ' nm'])
