%COLORMAPS3C lets you create your own colormaps.
% COLORMAPS3C is able to create a custom colormap by specifying the 3 colors
% that you want to use in your colormap. The first will be at the bottom of
% the scale, the second will be the intermediate and the third will be at the
% top of the scale.
% If you execute this program after that a figure with you own data has been
% generated, the COLORMAPS3C will apply automatically the new colormap
% generating a new colorbar.
% See also COLORBAR, COLORMAP.

%   Danilo Botta 17-05-04

sc=1;
echo off;


while sc==1,
    A=uisetcolor([0 0 0], 'Select first color:');
    X=uisetcolor([0.5 0.5 0.5], 'Select second color:');
    B=uisetcolor([1 1 1], 'Select third color:');
    

    C1=(X-A);
    C2=(B-X);

    for i=1:128,
        special(i,:)=[(A(1,1)+(C1(1,1)*i)/128) (A(1,2)+(C1(1,2)*i)/128) (A(1,3)+(C1(1,3)*i)/128)];
    end    

    for i=129:255,
        special(i,:)=[X(1,1)+(C2(1,1)*(i-128)/127) X(1,2)+(C2(1,2)*(i-128)/127) X(1,3)+(C2(1,3)*(i-128)/127)];
    end    

    colormap(special);
    colorbar;
    
    clear sc;
      sc=menu('Do you wish to repeat COLORMAPS3C?', 'Repeat', 'End');
end