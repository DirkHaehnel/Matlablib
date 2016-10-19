function y = crosscorr(x,xx) 
% CrossCorr(x) calculates the autocorrelation of a vector x

if size(x,1)==1 | size(x,2)==1
   x = x(:); len = length(x);
else
   len = size(x,1);
end
a = x;% - ones(len,1)*mean(x);
aa = xx;% - ones(len,1)*mean(x);
z = fft([a; zeros(size(a))]); 
zz = fft([aa; zeros(size(aa))]); 
z = real(ifft(zz.*conj(z)));
y = z(1:len,:)./((len:-1:1)'*ones(1,size(a,2)));

