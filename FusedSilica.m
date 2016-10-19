function [n,k] = FusedSilica(Wavelength)

if Wavelength < 400 
    display('Outside of accepted ranges')
    return
elseif Wavelength > 800
    display('Outside of accepted ranges')
    return
end

Wvlg = Wavelength/1000;

AZero =  2.104025406;
AOne  = -1.456000330e-04;
ATwo  = -9.049135390e-03;
AThree=  8.801830992e-03;
AFour =  8.435237228e-05;
AFive =  1.681656789e-06;
ASix  = -1.675425449e-08;
ASevn =  8.326602461e-10;

nSquare = AZero;
nSquare = nSquare + AOne * Wvlg.^4;
nSquare = nSquare + ATwo * Wvlg.^2;
nSquare = nSquare + AThree * Wvlg.^-2;
nSquare = nSquare + AFour * Wvlg.^-4;
nSquare = nSquare + AFive * Wvlg.^-6;
nSquare = nSquare + ASix * Wvlg.^-8;
nSquare = nSquare + ASevn * Wvlg.^10;

n = nSquare.^0.5;
k = zeros(size(n));
