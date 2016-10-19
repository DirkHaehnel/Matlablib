clear all;
close all;

% Manually choose files

[filename, pathname] = uigetfile('*.bh','Choose ASCII data file:');
%filename


% Measurement Parameters

framenbr = 1; % frame number
ADC = 256;
xdimension = 128;
ydimension = 128;

framesize = ADC*xdimension*ydimension;
framesizebytes = 2*ADC*xdimension*ydimension;
framestartpos = (framenbr-1)*framesizebytes+1;

fid = fopen(filename);
fseek(fid,framestartpos,'bof');
prerawdata = fread(fid,framesize,'uint16'); %use inf to verify
rawdata = prerawdata ./ 25.6; % 256? - discrepency btw matlab code, vistavision export and B&H SPCM export
clear prerawdata;
fclose(fid);

% Create image

imageintensity = zeros(xdimension,ydimension,'uint16');

for i = 1:(xdimension-1)
    for j = 1:(ydimension-1)
        cut = rawdata( ((i-1)*ydimension*ADC + (j-1)*ADC + 1) : ((i-1)*ydimension*ADC + j*ADC) ); % cut histogram for pixel (i,j)
        for k = 1:ADC
            imageintensity(i,j) = imageintensity(i,j) + cut(k);
        end
    end
end

h = imagesc(imageintensity);
colorbar;
colormap(gray);
clear imageintensity;

% Create histograms array of whole image

lifetimearray = zeros(xdimension,ydimension,ADC);
for i = 1:(xdimension-1)
   for j = 1:(ydimension-1)
       cut = rawdata( ((i-1)*ydimension*ADC + (j-1)*ADC + 1) : ((i-1)*ydimension*ADC + j*ADC) ); % cut histogram for pixel (i,j)
       lifetimearray(i,j,:) = cut; %3D array of histogram per image pixel position 
   end
end
%clear rawdata;


%%
% Lifetime evaluation

t = 1:size(lifetimearray,3); 
ind = t>50 & t<200; 
y = squeeze(sum(sum(lifetimearray,1),2)); 
p = Simplex('ExpFun',[10 20],[0 0],[],[],[],t(ind),y(ind),1);
[err, c, zz, z] = ExpFun(p,t(ind),y(ind),1,1);
clf
semilogy(t,y,t(ind),zz.*(ones(sum(ind),1)*c'))

for i = 1:(xdimension-1)
   for j = 1:(ydimension-1)
       amp(i,j,:) = lsqnonneg(zz,squeeze(lifetimearray(i,j,ind))); %fitting the lifetime patterns
   end
end

mim((amp(:,:,2)*p(1)+amp(:,:,3)*p(2))./sum(amp(:,:,2:3),3),sum(amp,3))
%%


return

% Access and output lifetime histogram of selected pixel

%xpos = 1; % x position of interest
%ypos = 1; % y position of interest
 disp('--->  select x0, y0')
 [X,Y] = ginput (1); %Get coordinates selected by mouse
 xpos = round(X);
 ypos = round(Y);
prehist = rawdata( ((xpos-1)*ydimension*ADC + (ypos-1)*ADC + 1) : ((xpos-1)*ydimension*ADC + ypos*ADC) ); % lifetimearray for point of interest
histm = squeeze(prehist);
for t = 1:ADC
    time(t) = 50/ADC*t;
end
figure; plot (time,histm);
save('histogram.txt','histm','-ASCII','-tabs');
% clear unnecessary;

% Fit Curve

% log scale y axis
% find max
% tail-fitting

%clear all;
    
        
        
        
    
