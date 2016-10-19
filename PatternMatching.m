% reading in the image
if exist('im')==0
    im = double(imread('c:\Joerg\Eudora\attach\molecule example 1.png'));
    im = im(100:end-101,100:end-100);
end

% defining the optical parameters of the imaging setup
z = 0;
NA = 1.40;
n0 = 1.52;
n = 1.49;
n1 = 1.0;
d0 = [];
d = 0.01;
d1 = [];
lamem = 0.57;
mag = 140;
focus = 0.65;
atf = [];
ring = [];
pixel = 16;
pic = 0;
be_res = [];
al_res = [];
nn = [];

% calculating the ideal images
model = PatternGeneration(z, NA, n0, n, n1, d0, d, d1, lamem, mag, focus, atf, ring, pixel, nn, be_res, al_res, pic);

% finding the patterns
bck = Disk((size(model.mask,1)-1)/2);
[err, bim, cim, sim, xc, yc, bc, cc, sc, len, imm] = FindPattern(im,model.mask,bck,bck);

% showing the results
CombineImages(cat(3,im,imm),1,2)


