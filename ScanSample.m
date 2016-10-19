function ScanSample(name, xoff, yoff, xc, yc, dt, nframes);

xc = xc + xoff;
yc = yc + yoff;

fprintf('Total numer of frames: %03d',nframes*size(xc')+1);

for i = 1:size(xc')
    fname = sprintf('%s_%03d',name, i);
    TH_ReadPixel(fname, xc(i), yc(i), dt, nframes);
end;

fname = sprintf('%s_%03d.t3r',name, 0);
TH_ReadPixel(fname, 0, 0, dt, nframes);

