function [err, c, z, cx, zx] = Riglerho(p, t, y, m, yx, delta, bld)

if nargin>4 && ~isempty(yx)
   if isempty(delta)
       delta = p(3);
   end
end
if length(t)<length(y)
    c = t;
    t = y(:);
    p = p(:)';
    if nargin<4 || isempty(m)
        m = 2;
    end
    switch m
        case 2
            z = 1./(p(1)^2+4*t)./sqrt(p(2)^2+4*t);
            if nargin>4 && ~isempty(yx)
                zx = z.*exp(-delta.^2./(p(1)^2+4*t));
            end
        case 3
            z = 1./(6*p(1)^4+64*p(1)^2*t+128*t.^2)./sqrt(6*p(2)^4+64*p(2)^2*t+128*t.^2);
            if nargin>4 && ~isempty(yx)
                zx = z.*exp(-4*delta.^2*(p(1)^2+4*t)./(3*p(1)^4+32*p(1)^2*t+64*r.^2));
            end
        case 4
            z = 1./(2*p(1)^6+40*p(1)^4*t+192*p(1)^2*t.^2+256*t.^3)./sqrt(2*p(2)^6+40*p(2)^4*t+192*p(2)^2*t.^2+256*t.^3);
            if nargin>4 && ~isempty(yx)
                zx = z.*exp(-delta.^2*(1./(2*p(1)^2+8*t) + (p(1)^2+4*t)./(p(1)^4+16*p(1)^2*t+32*t.^2)));
            end
    end
    err = z*c;
else
    t = t(:);
    y = y(:);
    p = p(:)';
    if nargin<4 || isempty(m)
        m = 2;
    end
    switch m
        case 2
            z = 1./(p(1)^2+4*t)./sqrt(p(2)^2+4*t);
            if nargin>4 && ~isempty(yx)
                zx = z.*exp(-delta.^2./(p(1)^2+4*t));
            end
        case 3
            z = 1./(6*p(1)^4+64*p(1)^2*t+128*t.^2)./sqrt(6*p(2)^4+64*p(2)^2*t+128*t.^2);
            if nargin>4 && ~isempty(yx)
                zx = z.*exp(-4*delta.^2*(p(1)^2+4*t)./(3*p(1)^4+32*p(1)^2*t+64*r.^2));
            end
        case 4
            z = 1./(2*p(1)^6+40*p(1)^4*t+192*p(1)^2*t.^2+256*t.^3)./sqrt(2*p(2)^6+40*p(2)^4*t+192*p(2)^2*t.^2+256*t.^3);
            if nargin>4 && ~isempty(yx)
                zx = z.*exp(-delta.^2*(1./(2*p(1)^2+8*t) + (p(1)^2+4*t)./(p(1)^4+16*p(1)^2*t+32*t.^2)));
            end
    end
    c = lsqnonneg(z,y);
    z = z*c;
    if nargin>4 && ~isempty(yx)
        cx = lsqnonneg(zx,yx);
        zx = cx*zx;
    end
    if nargin>6 && ~isempty(bld)
        if nargin>4 && ~isempty(yx)
            semilogx(t,y,'o',t,z,t,yx,'o',t,zx, 'markersize', 2.5); drawnow;            
        else
            semilogx(t,y,'o',t,z, 'markersize', 2.5); drawnow;
        end
    end
    err = sum(sum((y-z).^2./abs(z)));
    if nargin>4 && ~isempty(yx)
        err = err + sum((yx-zx).^2./abs(zx));
    end
end
