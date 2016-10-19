function [irf,y] = jw_simtcspc(A,T)
%The function jw_simtcspc simulates time-resolved data for testing TSCPCfit.
% It is called by: [irf,y] = jw_simtcspc(A,T).
% The function arguments are:
% A 	= 	Amplitudes A for the decay curve A*exp(-t/tau) of each species
% T 	= 	Lifetimes tau for the decay curve A*exp(-t/tau) of each species
%
%Internal parameters are:
% period p between laser pulses is p = 256 channels
% time duration dt of 1 channel is dt = 1 channel
%
%The return parameters are:
% irf = instrument response (mean = channel 30, fwhm = 10 channels)
% y = data WITHOUT ERROR IN THE PRESENT VERSION
% (c) 2007 John Williams

%Simulated instrument response according to page 375 in
% Enderlein and Erdmann (1997) Optics Communications 134:371
% "The instrument response function was modeled by a Gaussian curve witih
% full width at half maxiumu of 10, its maximum position in the 30th
% channel, and a peak height of 10^3 counts.
% THE MATH WAS WORKED OUT BY JOHN W IN THE MAPLE FILE: SimLifetime.mw

% time interval is length of 1 channel, there are 256 channels
dt = 1; 

% form the irf (K) for 1 cycle comprising 256 channels; see Maple
% SimLifetime.mw for the math
x = 1:256;
x4 = 1:4*256;
K = 1000 * exp(-1/25 * (x-30).^2 * log(2));
K4 = [K K K K];

kInit = 3*256 + 1; % init at channel 1 of the 4th (last) cycle
sig = []; % initialize output array 
for k = kInit:kInit+255 % for each channel k in the 4th (last) cycle 
    sj = 0;
    for j = 0:k-1 % integrate signal over each channel j going back 4 cycles
        si = 0;
        for i=1:length(A) % for each species i
            si = si + ( A(i)* exp(-j*dt/T(i)) * K4(k-j) );
        end
        sj = sj + si;
    end
    sig = [sig sj];
end
irf = K;
y = sig;
sig4 = padarray(sig,[0,3*256],'pre');
plot(x4,K4,'r',x4,sig4,'b')

        
            
    
    
