% Regenbogen
close

phi = 0:pi/500:2*pi;
n = 1.33;

thetav = 0:pi/180:pi/2;

res = zeros(2,2*length(thetav));
for j=1:length(thetav)
    theta = thetav(j);

    al = asin(sin(theta)/n);
    om1 = theta + pi - 2*al; % reflexion point
    om2 = om1 + pi - 2*al; % exit point

    bigtheta = om2 + theta;

    plot(cos(phi),sin(phi),'g','Linewidth',2);
    hold on
    %plot([-cos(theta) -cos(theta)*1.4],[sin(theta) sin(theta)*1.4],'b',[-cos(om2) -cos(om2)-0.4],[sin(om2) sin(om2)],'b')
    plot([-cos(theta)-1.5 -cos(theta)],sin(theta)*[1 1],'b')
    plot([-cos(theta) -cos(om1)], [sin(theta) sin(om1)],'b')
    plot([-cos(om1) -cos(om2)], [sin(om1) sin(om2)],'b')
    res(:,2*j-1:2*j) = [[-cos(om2) -cos(om2)-1.5*cos(bigtheta)]', [sin(om2) sin(om2)+1.5*sin(bigtheta)]'];
    for k=1:j-1
        plot(res(:,2*k-1),res(:,2*k))
    end
    plot(res(:,2*j-1),res(:,2*j),'b')
    hold off
    axis image
    axis off
    axis([-3 1.2 -2 1.2])
    drawnow
end