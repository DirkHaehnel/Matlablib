shift_size_xv=0.5:0.05:1.5;
shift_size_yv=-(0.5:0.05:1.5);

for jjx=1:length(shift_size_xv)
    for jjy=1:length(shift_size_yv)
        shift_size_x = shift_size_xv(jjx);
        shift_size_y = shift_size_yv(jjy);
        BerndFourier;
        [mx0(jjx,jjy),my0(jjx,jjy),wx0(jjx,jjy),wy0(jjx,jjy)] = Gauss2D(im0(a,b));
        [mx1(jjx,jjy),my1(jjx,jjy),wx1(jjx,jjy),wy1(jjx,jjy)] = Gauss2D(im(a,b));
        %[mx2(jjx,jjy),my2(jjx,jjy),wx2(jjx,jjy),wy2(jjx,jjy)] = Gauss2D(res(a,b));
    end
end
