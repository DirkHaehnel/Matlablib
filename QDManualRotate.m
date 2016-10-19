load('D:\Doc\Patra\SauerDots\picks.mat')

for j=1:size(imm,3)
    surf(imm(:,:,j)); view([0 90]); axis image; shading flat; colormap hot; cameratoolbar('setmode','roll')
    pause
    a = view;
    phi(j) = atan(a(1,2)/a(1,1))/pi*180;
end

