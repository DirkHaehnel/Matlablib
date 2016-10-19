function y = Autocorr(x) 
% Autocorr(x) calculates the autocorrelation of a vector x

if size(x,1)==1 || size(x,2)==1
   x = x(:); len = length(x);
else
   len = size(x,1);
end
a = x - ones(len,1)*mean(x);
z = fft([a; zeros(size(a))]); 
z = real(ifft(abs(z).^2));
y = z(1:len,:)./((len:-1:1)'*ones(1,size(a,2)));

