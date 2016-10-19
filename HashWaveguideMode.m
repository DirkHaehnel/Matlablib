function [qp, qs, fp, fs] = HashWaveguideMode(n,d)

% function determines the in-plane wave vector components for the waveguide modes
% (c) Joerg Enderlein (2005)

n = n(:).'; d = d(:).';
dq = 1e-6; % determines the finesse of searching; if waveguide modes are closer than this spacing one has a problem 
qq = max([n(1) n(end)])+dq/2:dq:max(n)-dq/2;
d = [0 d 0];
j = length(n);
w1 = sqrt(n(j)^2-qq.^2);
w1(imag(w1)<0) = conj(w1(imag(w1)<0));
w0 = sqrt(n(j-1)^2-qq.^2);
w0(imag(w0)<0) = conj(w0(imag(w0)<0));
M11 = (w1./w0*n(j-1)/n(j)+n(j)/n(j-1))/2;
M12 = (-w1./w0*n(j-1)/n(j)+n(j)/n(j-1))/2;
M21 = M12;
M22 = M11;
for j=length(n)-1:-1:2
    w1 = sqrt(n(j)^2-qq.^2);
    w1(imag(w1)<0) = conj(w1(imag(w1)<0));
    w0 = sqrt(n(j-1)^2-qq.^2);
    w0(imag(w0)<0) = conj(w0(imag(w0)<0));
    M11 = exp(-i*w1*d(j)).*M11;
    M12 = exp(-i*w1*d(j)).*M12;
    M21 = exp(i*w1*d(j)).*M21;        
    M22 = exp(i*w1*d(j)).*M22;                
    N11 = (w1./w0*n(j-1)/n(j)+n(j)/n(j-1))/2;
    N12 = (-w1./w0*n(j-1)/n(j)+n(j)/n(j-1))/2;
    N21 = N12;
    N22 = N11;
    tmp11 = N11.*M11+N12.*M21;
    tmp12 = N11.*M12+N12.*M22;
    M21 = N21.*M11+N22.*M21;
    M22 = N21.*M12+N22.*M22;
    M11 = tmp11;
    M12 = tmp12;
end
fp = M11;

j = length(n);
w1 = sqrt(n(j)^2-qq.^2);
w1(imag(w1)<0) = conj(w1(imag(w1)<0));
w0 = sqrt(n(j-1)^2-qq.^2);
w0(imag(w0)<0) = conj(w0(imag(w0)<0));
M11 = (w1./w0+1)/2;
M12 = (-w1./w0+1)/2;
M21 = M12;
M22 = M11;
for j=length(n)-1:-1:2
    w1 = sqrt(n(j)^2-qq.^2);
    w1(imag(w1)<0) = conj(w1(imag(w1)<0));
    w0 = sqrt(n(j-1)^2-qq.^2);
    w0(imag(w0)<0) = conj(w0(imag(w0)<0));
    M11 = exp(-i*w1*d(j)).*M11;
    M12 = exp(-i*w1*d(j)).*M12;        
    M21 = exp(i*w1*d(j)).*M21;        
    M22 = exp(i*w1*d(j)).*M22;                
    N11 = (w1./w0+1)/2;
    N12 = (-w1./w0+1)/2;
    N21 = N12;
    N22 = N11;
    tmp11 = N11.*M11+N12.*M21;
    tmp12 = N11.*M12+N12.*M22;
    M21 = N21.*M11+N22.*M21;
    M22 = N21.*M12+N22.*M22;
    M11 = tmp11;
    M12 = tmp12;
end
fs = M11;

t = 1:length(fp);
t = t(fp(1:end-1).*fp(2:end)<0);
qp = [];
for j=1:length(t)
    ind = max([1 t(j)-2]):min([length(qq) t(j)+3]);
    qp(j) = interp1(fp(ind), qq(ind), 0, 'spline');
end

t = 1:length(fs);
t = t(fs(1:end-1).*fs(2:end)<0);
qs = [];
for j=1:length(t)
    ind = max([1 t(j)-2]):min([length(qq) t(j)+3]);    
    qs(j) = interp1(fs(ind), qq(ind), 0, 'spline');
end

