close all
scrsz = get(0,'ScreenSize'); 
n = 7;

figure('Name',['Gaussian fits x ' fname], 'NumberTitle','off');
plot(t,hx,'.',t,htx);
figure('Name',['Gaussian fits y ' fname], 'NumberTitle','off');
plot(t,hy,'.',t,hty);
figure('Name',['Gaussian fits z ' fname], 'NumberTitle','off');
plot(t,hz,'.',t,htz);

figure('Name',['Parameter fits ' fname], 'NumberTitle','off', 'Position',[scrsz(3)/2-500 scrsz(4)/2-300 1000 400]);
subplot(131);
[err, cx, gausstx] = OneMExpFun(kx, 2.^(0:n-1), gaussx.^2);
xlabel('Time');
title('x-coordinate');
ylabel('Width of Gaussian');
ax = axis; ax(4) = ax(4)-ax(3);
text(0.7*2^(n-1),ax(3) + 0.2*ax(4),['kx = ' num2str(kx)]);
text(0.7*2^(n-1),ax(3) + 0.1*ax(4),['cx = ' num2str(cx(end))]);
subplot(132);
[err, cy, gaussty] = OneMExpFun(ky, 2.^(0:n-1), gaussy.^2);
xlabel('Time');
title('y-coordinate');
ax = axis; ax(4) = ax(4)-ax(3);
text(0.7*2^(n-1),ax(3) + 0.2*ax(4),['ky = ' num2str(ky)]);
text(0.7*2^(n-1),ax(3) + 0.1*ax(4),['cy = ' num2str(cy(end))]);
subplot(133);
[err, cz, gausstz] = OneMExpFun(kz, 2.^(0:n-1), gaussz.^2);
xlabel('Time');
title('z-coordinate');
ax = axis; ax(4) = ax(4)-ax(3);
text(0.7*2^(n-1),ax(3) + 0.2*ax(4),['kz = ' num2str(kz)]);
text(0.7*2^(n-1),ax(3) + 0.1*ax(4),['cz = ' num2str(cz(end))]);

