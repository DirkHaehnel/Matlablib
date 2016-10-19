% Defocus imaging
% calculation under the equation in JCPA by Dr J. Enderlein
%

clear;

%-----------------------------------------------------
% experimental conditions
%-----------------------------------------------------
%
% best for demo: NA=1.2, dz=1.4, z=0.08, lambda=0.6
%
%



% defocusing shift of objective toward the sample (in um)
dz = -0.9;

% location of dye blow the surface (in um)
z = 1/1000* 40;

% numerical aperture of objective
NA = 1.35;

% magnification of microscope
mag = 100;

% refractive index of layer above the sample
n0 = 1.0;

% refractive index of the sample and glass
n1 = 1.52;

	%   n0
	% ----------
	%       } z
	%   n1  *

% pixel size of CCD chip in um and (pixel number-1)/2
PixelSize = 6.45;
PixelRes  = 20;

% number of division in calculation of integral
ND = 750;

% wavelength and wavenumber of the probe light (in um)
lambda = 0.6;
%lambda = repetition;
k = 2*pi/lambda;

% polar angle of the dye molecule (in degree)
q = 90;
Theta = q*pi/180;



%------------------------------------------------------------------
% Fresnel reflection and transmission coefficients
% 	Rp, Rs: reflection for p- and s-polarization
% 	Tp, Ts: transmission for p- and s-polarization
% 	n1,n2: refractive indeces: direction of 1=>2


function [Rp, Rs, Tp, Ts] = Fresnel(theta, n1, n2);
	if (theta==0.)
		Rp = (-n1+n2)./(n1+n2);
		Rs =  (n1-n2)./(n1+n2);
		Tp = 2.*n1./(n1+n2);
		Ts = 2.*n1./(n1+n2);
	elseif (theta == pi/2)
		Rp = 1;
		Rs = 1;
		Tp = 0;
		Ts = 0;
	else
		phi1 = theta;
		phi2 = asin(n1.*sin(phi1)./n2);
		Rp =  tan(phi1 - phi2)./tan(phi1 + phi2);
		Rs = -sin(phi1 - phi2)./sin(phi1 + phi2);
		Tp = 2.*cos(phi1).*sin(phi2)./sin(phi1 + phi2)./cos(phi1 - phi2);
		Ts = 2.*cos(phi1).*sin(phi2)./sin(phi1 + phi2);
	endif
endfunction

function [Rp, Rs, Tp, Ts] = Fresnelw(w1, n1, n2);

	w2 = sqrt(n2^2 - n1^2 + w1.^2);
	%w2(imag(w2)<0) = conj(imag(w2)<0);
	if imag(w2)<0
		w2 = conj(w2);
	endif
	
	Rp = (w1*n2^2 - w2*n1^2)./(w1*n2^2 + w2*n1^2);
	Rs = (w1 - w2)./(w1 + w2);
	Tp = 2*n1*n2*w1./(w1*n2^2 + w2*n1^2);
	Ts = 2*w1./(w1 + w2);
	
endfunction





%------------------------------------------------------------------
% p- and s-components of electric vectors for vertically (perp) 
% and horizontally (para) oriented dipoles.
%
function [Eperp_p, Epara_p, Epara_s] = EVector(sin_eta, k, z, n0, n1);

	cos_eta = sqrt(1 - sin_eta.^2);
	[Rp,Rs,Tp,Ts] = Fresnel(asin(sin_eta), n1, n0);
	
	% temporary (because the above function doesnt work well for theta==0)
		Rp(1) = (-n1+n0)./(n1+n0);
		Rs(1) =  (n1-n0)./(n1+n0);
		Tp(1) = 2.*n1./(n1+n0);
		Ts(1) = 2.*n0./(n1+n0);
	Eperp_p = sin_eta.*(exp(-i*n1*k*z*cos_eta) + Rp.*exp(i*n1*k*z*cos_eta));
	Epara_p = cos_eta.*(exp(-i*n1*k*z*cos_eta) - Rp.*exp(i*n1*k*z*cos_eta));
	Epara_s = cos_eta.*(exp(-i*n1*k*z*cos_eta) + Rs.*exp(i*n1*k*z*cos_eta));

endfunction


%-----------------------------------------------------
% main routine
%-----------------------------------------------------

% variables
% 	Eta0:  polar angle in the sample space
% 	Eta1:  polar angle in the image space
% 	dEta1: difference of Eta1
% 	Eta1_max: highest angle of Eta1

[x,y] = meshgrid(-PixelRes:PixelRes, -PixelRes:PixelRes);
 rho  = sqrt(x.^2 + y.^2)*PixelSize;
 psi  = arg(x+i*y);
 img  = zeros(length(x));
 img_para = zeros(length(x));
 img_perp = zeros(length(x));
 img_cros = zeros(length(x));

Eta1_max = asin(NA/mag);
dEta1 = Eta1_max/ND;
Eta1 = (0: dEta1: Eta1_max);
sin_Eta1  = sin(Eta1);
cos_Eta1  = cos(Eta1);

% Eta0 = asin(mag.*sin(Eta1)/n1);
% cos_Eta0  = cos(Eta0);
% sin_Eta0  = sin(Eta0);

sin_Eta0  = 1/n1*mag*sin_Eta1;
cos_Eta0  = sqrt(1 - sin_Eta0.^2);


col_rho = rho(:);
col_psi = psi(:);
col_img = img(:);
col_img_para = img_para(:);
col_img_perp = img_perp(:);
col_img_cros = img_cros(:);

length_img = length(col_img);

[Eperp_p, Epara_p, Epara_s] = EVector(sin_Eta0, k, z, n0, n1);

for ind = 1:length_img
	cos_psi  = cos(col_psi(ind));
	sin_psi  = sin(col_psi(ind));
	cos_2psi = cos(2*col_psi(ind));
	sin_2psi = sin(2*col_psi(ind));

	besselarg = k*col_rho(ind).*sin_Eta1;
	J0 = besselj(0, besselarg);
	J1 = besselj(1, besselarg);
	J2 = besselj(2, besselarg);
	
	e_x_para = i/2.*(cos_Eta1.*(J0-J2*cos_2psi).*Epara_p + (J0+J2*cos_2psi).*Epara_s);
	e_y_para = i/2.*(-cos_Eta1.*J2.*sin_2psi.*Epara_p + J2.*sin_2psi.*Epara_s);
	b_x_para = i/2.*(-cos_Eta1.*J2.*sin_2psi.*Epara_s + J2.*sin_2psi.*Epara_p);
	b_y_para = i/2.*(cos_Eta1.*(J0+J2*cos_2psi).*Epara_s + (J0-J2*cos_2psi).*Epara_p);
	
	e_x_perp = i.*cos_Eta1.*J1.*Eperp_p.*cos_psi;
	e_y_perp = i.*cos_Eta1.*J1.*Eperp_p.*sin_psi;
	b_x_perp = i.*J1.*Eperp_p.*(-sin_psi);
	b_y_perp = i.*J1.*Eperp_p.*cos_psi;
	
	Ex_para = (sin_Eta1.*sqrt(cos_Eta1./n1./cos_Eta0).*e_x_para.*exp(i*k*dz*cos_Eta0)) *ones(length(Eta1),1) *dEta1;
	Ey_para = (sin_Eta1.*sqrt(cos_Eta1./n1./cos_Eta0).*e_y_para.*exp(i*k*dz*cos_Eta0)) *ones(length(Eta1),1) *dEta1;
	Bx_para = (sin_Eta1.*sqrt(cos_Eta1./n1./cos_Eta0).*b_x_para.*exp(i*k*dz*cos_Eta0)) *ones(length(Eta1),1) *dEta1;
	By_para = (sin_Eta1.*sqrt(cos_Eta1./n1./cos_Eta0).*b_y_para.*exp(i*k*dz*cos_Eta0)) *ones(length(Eta1),1) *dEta1;

	Ex_perp = (sin_Eta1.*sqrt(cos_Eta1./n1./cos_Eta0).*e_x_perp.*exp(i*k*dz*cos_Eta0)) *ones(length(Eta1),1) *dEta1;
	Ey_perp = (sin_Eta1.*sqrt(cos_Eta1./n1./cos_Eta0).*e_y_perp.*exp(i*k*dz*cos_Eta0)) *ones(length(Eta1),1) *dEta1;
	Bx_perp = (sin_Eta1.*sqrt(cos_Eta1./n1./cos_Eta0).*b_x_perp.*exp(i*k*dz*cos_Eta0)) *ones(length(Eta1),1) *dEta1;
	By_perp = (sin_Eta1.*sqrt(cos_Eta1./n1./cos_Eta0).*b_y_perp.*exp(i*k*dz*cos_Eta0)) *ones(length(Eta1),1) *dEta1;
	
	col_img_para(ind) = Ex_para*conj(By_para) - conj(Bx_para)*Ey_para;
	col_img_perp(ind) = Ex_perp*conj(By_perp) - conj(Bx_perp)*Ey_perp;
	col_img_cros(ind) = Ex_perp*conj(By_para) + Ex_para*conj(By_perp) - (Ey_perp*conj(Bx_para) + Ey_para*conj(Bx_perp));
	
endfor

img_para(:) = img_para(:) + col_img_para;
img_perp(:) = img_perp(:) + col_img_perp;
img_cros(:) = img_cros(:) + col_img_cros;
img = real(img_para.*(sin(Theta).^2) +img_perp.*(cos(Theta).^2) +img_cros.*sin(Theta).*cos(Theta));


%surf(img);
%plot(x(1,:),img_para(PixelRes+1,:), x(1,:),img_perp(PixelRes+1,:), x(1,:),img_cros(PixelRes+1,:));
%plot(y(:,1),img_para(:,PixelRes+1), y(:,1),img_perp(:,PixelRes+1), y(:,1),img_cros(:,PixelRes+1));
plot(x(1,:),img(PixelRes+1,:), y(:,1),img(:,PixelRes+1));


