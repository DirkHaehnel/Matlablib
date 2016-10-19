% RGB VALUES FOR VISIBLE WAVELENGTHS
%
% Adapted from fortran to MatLab by B. Müller, 
% using approximation from Dan Bruton (astro@tamu.edu)
% "Approximate RGB values for Visible Wavelengths, 1996"
% http://www.physics.sfasu.edu/astro/color.html
%
% Rainbow: clf; for j=1:251 patch([398.5+j 398.5+j 399.5+j 399.5+j], 20*[0 1 1 0],nm2rgb(400+j-1),'edgecolor','none'); end; axis image; axis off

function [ rgb ] = nm2rgb( nm )
%NM2RGB Summary of this function goes here
%   Detailed explanation goes here

if (nm >= 380) && (nm <= 780)

    if (nm >= 380) && (nm <= 440)
        r = -1 * (nm-440)/(440-380);
        g = 0;
        b = 1;
    end;
 
    if (nm > 440) && (nm <= 490)
        r = 0;
        g = (nm-440)/(490-440);
        b = 1;
    end;
 
    if (nm > 490) && (nm <= 510)
        r = 0;
        g = 1;
        b = -1*(nm-510)/(510-490);
    end;

    if (nm > 510) && (nm <= 580)
        r = (nm-510)/(580-510);
        g = 1;
        b = 0;
    end;

    if (nm > 580) && (nm <= 645)
        r = 1;
        g = -1*(nm-645)/(645-580);
        b = 0;
    end;

    if (nm > 645) && (nm <= 780)
        r = 1;
        g = 0;
        b = 0;
    end;
        
else
    r = 0;
    g = 0;
    b = 0;
end;
if  nm >= 700
  SSS = 0.3+0.7*(780-nm)/(780-700);
else
  if nm <= 420
     SSS = 0.3+0.7*(nm-380)/(420-380);
  else
     SSS = 1;
  end
end

rgb = [r*SSS g*SSS b*SSS];
end
