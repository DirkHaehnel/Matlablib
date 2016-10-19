shift_size_vec = (-0.5:0.1:0.5)*pixel_size/step_size;

for k=1:length(shift_size_vec)
    shift_size = shift_size_vec(k);
    im = zeros(nx,ny);
    for jx=1:sx
        for jy = 1:sy
            im = im + interp2(x,y,squeeze(imraw(jx,jy,:,:)),x+xx(jx)*shift_size,y+yy(jy)*shift_size,'cubic',0);
        end
    end

    mim(cat(3,squeeze(sum(sum(imraw,1),2)),im));
    title(k)
    eval(['print -dpng -r300 BerndTestPP' mint2str(k,2)])
end

