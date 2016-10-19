m = 2;

mask = diag(ones(2*m,1),-1)/2;
mask = mask + mask' + diag(ones(2*m+1,1));

tst = mconv2(tmp,mask);
tst = max(tst,mconv2(tmp,fliplr(mask)));
mask = zeros(2,2*m+1); 
mask(1,:) = 1;
mask(2,:) = 1/2;

tst = max(tst,mconv2(tmp,mask));
tst = max(tst,mconv2(tmp',mask)');
mim(tmp,tst);