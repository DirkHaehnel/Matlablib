% Gustafsson figure


x = -10:0.01:10; 
y0 = exp(-2*x.^2/0.2^2);
t = abs(x)<1;

sat = [1e-3 10.^(-1:0.1:1)];
ss = [0 10.^(-2:0.01:1)];

if 0
    xx = 2*pi/10*100*x;
    tt = abs(xx)<150;
    for j = 1:length(sat)
        y = 1 - exp(-sat(j)*y0); y = y/max(y);

        subplot(2,3,1)
        plot(ss, 1 - exp(-ss), sat(j), 1 - exp(-sat(j)), 'o')
        xlabel('excitation intensity'); ylabel('fluorescence intensity')

        subplot(2,3,2);
        plot(x(t),y(t));
        xlabel('position [\mum]'); ylabel('fluorescence intensity')

        yy = abs(fftshift(fft(y))); yy = yy/max(yy);
        subplot(2,3,3);
        plot(xx(tt),yy(tt));
        axis tight
        xlabel('frequency [1/\mum]'); ylabel('power')

        eval(['print -dpng -r300 tmp' mint2str(j,2)])
        eval(['print -dpng -r300 tmp' mint2str(2*length(sat)+1-j,2)])
    end
end


if 0
    close all
    xx = 0:0.5:50; 
    tau = 5;
    sat = 4;
    for j = 1:length(xx)
        y = sat*y0.^2./(1+sat*y0).*exp(-xx(j)*y0/tau) + y0./(1+sat*y0);

        plot(x(t),y(t));
        hold on
        plot(x(t),y0(t),'b',x(t),y0(t)./(1+sat*y0(t)),'c','linewidth',2);
        hold off
        xlabel('position [\mum]'); ylabel('fluorescence intensity')
        axis([-0.5 0.5 0 1])
        text(0.1,0.8,['time = ' mnum2str(xx(j),2,1) ' \mus'],'units','normal')

        drawnow
        eval(['print -dpng -r300 tmp' mint2str(j,2)])
    end
end