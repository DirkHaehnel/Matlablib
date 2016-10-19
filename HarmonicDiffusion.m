% Pogram for modling diffusion of a particle within a harmonic potential 

N = 50;
k = 1;
D = 10;

% initial conditions:
x = -(N/2-0.5):N/2;
%v = [zeros(N/2); diag(ones(1,N/2))];

v = zeros(N,1);
v(N/2) = 1;

one = ones(N,1);
D = expm(D*(diag(one(1:N-1),-1) + diag(one(1:N-1),1) -2*diag(one)));

t = 1;
for j=1:100
   v = D*v;
   v(1,:) = 0;
   v(end,:) = 0;
   t = t+1;
   plot(x/sqrt(k*t),v*sqrt(k*t)); axis([min(x/sqrt(k*t)) max(x/sqrt(k*t)) 0 1]); drawnow
end

   