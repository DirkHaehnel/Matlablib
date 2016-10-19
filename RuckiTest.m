% program RuckiTest: Kann man einen evanescenten Kollektor aus spherischen
% Flächen basteln?

R = 1e4;
theta = (asin(1.33/1.51):pi/1e2:pi/2)';

Lv = 1.02*R:1e1:1.04*R;

n0 = [sin(theta) 0*theta -cos(theta)];
r0 = 0*n0;
    
[r1, n1] = Reflection(r0, n0, [], [], R, -R*sin(asin(1.45/1.51)));

for j=1:length(Lv)
    [r2, n2] = Refraction(r1, n1, [], [], Lv(j), -1 + r1(1,3) - Lv(j)*cos(asin(r1(1,1)/Lv(j))), 1.51);
    n2(abs(imag(r2(:,1)))>0,:) = [];
    r2(abs(imag(r2(:,1)))>0,:) = [];
    lam = -r2(:,1)./n2(:,1);
    r3 = r2 + (lam*[1 1 1]).*n2;
    plot3(r1(:,1),r1(:,2),r1(:,3),r2(:,1),r2(:,2),r2(:,3),r3(:,1),r3(:,2),r3(:,3)); axis image; pause
    %plot([r0(:,1) r1(:,1)]',[r0(:,3) r1(:,3)]'); axis image
    %plot(theta,acos(-sum(direx.*n1,2))/pi*180)
    res(j) = std(r3(:,3));
    
    % that's not a good criterion; one should ask: how much spherical is the
    % wavefront 
end
