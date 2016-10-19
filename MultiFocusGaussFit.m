function [wx, wy, mx, my, om] = MultiFocusGaussFit(tim);

for j=1:3
    for k=1:2
        subplot(2,3,(k-1)*3+j)
        [mx(j,k) my(j,k) wx(j,k) wy(j,k) om(j,k)] = GaussEllipse([], [], tim(:,:,j,k));
    end
end

figure
phi = 0:pi/100:2*pi;
plot(mx,my,'s');
hold on
for j=1:3
    plot(mx(j,1)+wx(j,1)*cos(phi)*cos(om(j,1))+wy(j,1)*sin(phi)*sin(om(j,1)),my(j,1)-wx(j,1)*cos(phi)*sin(om(j,1))+wy(j,1)*sin(phi)*cos(om(j,1)),'r');
    plot(mx(j,2)+wx(j,2)*cos(phi)*cos(om(j,2))+wy(j,2)*sin(phi)*sin(om(j,2)),my(j,2)-wx(j,2)*cos(phi)*sin(om(j,2))+wy(j,2)*sin(phi)*cos(om(j,2)),'b');
end
axis image; axis([1 size(tim,2) 1 size(tim,1)]); 
hold off

