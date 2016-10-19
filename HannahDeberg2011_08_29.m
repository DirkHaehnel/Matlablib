z = 0; 
NA = 1.65; 
n0 = 1.78; 
n = 1.42; 
n1 = n; 
d0 = []; 
d = 0.0; 
d1 = []; 
lamem = 0.65; 
mag = 190; 
focus = 0; 
atf = []; 
ring = []; 
pixel = 16; 
pic = 0; 
be_res = [];
al_res = []; 
nn = 6;

model = PatternGeneration(z, NA, n0, n, n1, d0, d, d1, lamem, mag, focus, atf, ring, pixel, nn, be_res, al_res, pic);

% CombineImages(model.mask,15,15)

if ~(exist('im')==1)
    im = double(imread('c:\Joerg\Doc\Microscopy\Selvin\Hannah Deberg\CF633_ProlongGold_1.65NA.tif'));
    a = [145 254];
    b = [149 255];
    im = im(a(1):a(2),b(1):b(2));
end

bck = Disk((size(model.mask,1)-1)/2); 
sze = 1;
tsh = 1;
[err, bim, cim, sim, xc, yc, bc, cc, sc, len, imm] = FindPattern(im,model.mask,bck,bck,sze,tsh);

CombineImages(cat(3,im,imm),1,2); hold on; plot(xc,yc,'oc'); hold off