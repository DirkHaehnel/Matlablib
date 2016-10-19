function y = HellPhaseMod(x,om,nn)

if nargin<3
    nn = 10;
end
A = zeros(nn);
for j=1:nn
    switch j
        case 1,
            A(1,1) = (1+x);
            A(1,2) = x/2;
        case 2
            A(2,1) = x.*(1+x)/2./om.^2;
            A(2,2) = 1 + (8+x.*(16+9*x))/8./om.^2;
            A(2,3) = 3*x.*(1+x)/4./om.^2;
            A(2,4) = x.^2/8./om.^2;
        otherwise
            A(j,j-2) = 1/4/(j-2)/(j-1)*x.^2./om.^2;
            A(j,j-1) = (2*j-3)/(j-2)/(j-1)^2*x.*(1+x)./om.^2;
            A(j,j) = 1 + 1/2/((j-1)^2-1)*(2+x*(4+3*x))./om.^2-1/(j-1)^2/((j-1)^2-1)*(1+x).^2./om.^2;
            if j<nn
                A(j,j+1) = (2*j-1)/j/(j-1)^2*x.*(1+x)./om.^2;
            end
            if j<nn-1 
                A(j,j+2) = 1/4/j/(j-1)*x.^2./om.^2;
            end
    end
end
y = A\[1;zeros(nn-1,1)];

