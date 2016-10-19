function [tw  tmax wavelength dv] = CavityTransmission(met, nwater, d0, d1, dv, NA)

load metals
nglass = 1.52;
if NA>0
    theta = (0.5:ceil(asin(NA/1.33)/pi*180))/ceil(asin(NA/1.33)/pi*180)*asin(NA/1.33);
end 

ind = wavelength>=450;
wavelength = wavelength(ind);
eval([met '=' met '(ind);']);

% transmission calculation
tw = zeros(numel(wavelength),numel(dv));
for j=1:length(wavelength)
    eval(['nag = ' met '(wavelength==wavelength(j));'])
    for k=1:length(dv)
        nrm = 0;
        if NA>0
            for m=1:length(theta)
                [~, ~, tmp] = Fresnel(sqrt(nglass^2-1.33^2*sin(theta(m))^2),[nglass nag nwater nag nglass],[d0 dv(k) d1]/wavelength(j)*2*pi);
                tw(j,k) = tw(j,k) + sin(theta(m))*tmp;
                nrm = nrm + sin(theta(m));
            end
            tw(j,k) = tw(j,k)/nrm;
        else
            [~, ~, tw(j,k)] = Fresnel(nglass,[nglass nag nwater nag nglass],[d0 dv(k) d1]/wavelength(j)*2*pi);
        end
    end
end
tw = abs(tw).^2;
mpcolor(dv,wavelength,tw)
xlabel('cavity height (nm)')
ylabel('wavelength (nm)')

tmax = zeros(size(dv));
for k=1:length(dv)
    tmax(k) = wavelength(tw(:,k)==max(tw(:,k)));
end

