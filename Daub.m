function y=daub(x,isign);

% The function daub performs a Daubechie wavelet transform

[m,n]=size(x);

if (m>1 & n>1 & ~(rem(log2(m),1)==0)) | (m==1 & ~(rem(log2(n),1)==0)) | (n==1 & ~(rem(log2(m),1)==0))
   error('Length of vecor is not a power of 2');
end
if m==1 | n==1 x=x(:); end
n = size(x,1);

c = [(1+sqrt(3))/4/sqrt(2) (3+sqrt(3))/4/sqrt(2) (3-sqrt(3))/4/sqrt(2) (1-sqrt(3))/4/sqrt(2)];


cs = c;
cd = [c(4) -c(3) c(2) -c(1)];

tmp = spdiags(ones(n,1)*cs,0:3,n,n);
DWT = spdiags(ones(n,1)*cd,0:3,n,n);

DWT = [tmp(1:2:end,:); DWT(1:2:end,:)];
if isign>0
   for m=log2(n):-1:2
      tmp = x;
      x(1:2^m,:) = DWT([1:2^(m-1) (n/2+1):(n/2+2^(m-1))],1:2^m)*x(1:2^m,:);
      x(2^(m-1),:) = x(2^(m-1),:) + c(3:4)*tmp(1:2,:);
      x(2^m,:) = x(2^m,:) + [c(2) -c(1)]*tmp(1:2,:);
   end
else
   DWT = DWT';
   for m=2:log2(n)
      tmp = x;
      x(1:2^m,:) = DWT(1:2^m, [1:2^(m-1) (n/2+1):(n/2+2^(m-1))])*x(1:2^m,:);
      x(1,:) = x(1,:) + c([3 2])*tmp([2^(m-1) 2^m],:);
     	x(2,:) = x(2,:) + [c(4) -c(1)]*tmp([2^(m-1) 2^m],:);
   end
end
y=x;

bld = [y(1)*ones(1,n/4) y(2)*ones(1,n/4)];
tmp = zeros(1,n/2);
for m=2:log2(n)
   for j=1:2^(m-1)
      tmp(j:2^(m-1):n/2) = y(2^(m-1)+j);
   end
   bld = [bld; tmp];
end
figure; bld=pcolor(bld), set(bld,'edgecolor','none');