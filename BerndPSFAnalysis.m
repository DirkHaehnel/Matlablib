close all

fac = step_size*1e3;

mim(im0)
[a,b] = PickPixel;

if length(a)==1
    
    x = -17:17;
    py0 = Simplex('Gauss',[0 2],[],[],[],[],x,im0(a+x,b)/max(im0(a+x,b)),0,[],1,im0(a+x,b)/max(im0(a+x,b)));
    py1 = Simplex('Gauss',[0 2],[],[],[],[],x,im(a+x,b)/max(im(a+x,b)),0,[],1,im(a+x,b)/max(im(a+x,b)));
    py2 = Simplex('Gauss',[0 2],[],[],[],[],x,res(a+x,b)/max(res(a+x,b)),0,[],1,res(a+x,b)/max(res(a+x,b)));
    
    xi = 1.1*x(1):0.1:1.1*x(end);
    y0 = exp(-xi.^2/2/py0(2)^2);
    y1 = exp(-xi.^2/2/py1(2)^2);
    y2 = exp(-xi.^2/2/py2(2)^2);
    
    xi = fac*xi;
    plot(fac*(x-py0(1)),im0(a+x,b)/max(im0(a+x,b)),'o',xi,y0,'r',fac*(x-py1(1)),im(a+x,b)/max(im(a+x,b)),'vg',xi,y1,'g',fac*(x-py2(1)),res(a+x,b)/max(res(a+x,b)),'db',xi,y2,'b')
    axis([xi(1)/1.1 xi(end)/1.1 0 1.1])
    hold on; plot(xi,1/exp(2),':k'); hold off
    xlabel('position (nm)')
    ylabel('rel. intensity')
    
    figure
    
    px0 = Simplex('Gauss',[0 2],[],[],[],[],x,im0(a,b+x)/max(im0(a,b+x)),0,[],1,im0(a,b+x)/max(im0(a,b+x)));
    px1 = Simplex('Gauss',[0 2],[],[],[],[],x,im(a,b+x)/max(im(a,b+x)),0,[],1,im(a,b+x)/max(im(a,b+x)));
    px2 = Simplex('Gauss',[0 2],[],[],[],[],x,res(a,b+x)/max(res(a,b+x)),0,[],1,res(a,b+x)/max(res(a,b+x)));
    
    xi = 1.1*x(1):0.1:1.1*x(end);
    y0 = exp(-xi.^2/2/px0(2)^2);
    y1 = exp(-xi.^2/2/px1(2)^2);
    y2 = exp(-xi.^2/2/px2(2)^2);
    
    xi = fac*xi;
    plot(fac*(x-px0(1)),im0(a,b+x)/max(im0(a,b+x)),'o',xi,y0,'r',fac*(x-px1(1)),im(a,b+x)/max(im(a,b+x)),'vg',xi,y1,'g',fac*(x-px2(1)),res(a,b+x)/max(res(a,b+x)),'db',xi,y2,'b')
    axis([xi(1)/1.1 xi(end)/1.1 0 1.1])
    hold on; plot(xi,1/exp(2),':k'); hold off
    xlabel('position (nm)')
    ylabel('rel. intensity')
    
    [px0(2)/px2(2) py0(2)/py2(2)]
    
else
    
    a = a(1):a(2);
    b = b(1):b(2);
    
    [mx0,my0,wx0,wy0] = Gauss2D(im0(a,b));
    [mx1,my1,wx1,wy1] = Gauss2D(im(a,b));
    [mx2,my2,wx2,wy2] = Gauss2D(res(a,b));
    
    [wx0/wx1 wx1/wx2 wx0/wx2]
    [wy0/wy1 wy1/wy2 wy0/wy2]
    
end

