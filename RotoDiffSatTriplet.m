function y = RotoDiffSatTriplet(x,drot,kisc,kph)

nn = 10;
A = zeros(nn);
for j=1:nn
    if j>1 
        A(j,j-1) = (drot*(2*j-2)*(2*j-1)+kisc+kph)/(1+drot*(2*j-2)*(2*j-1)+kisc)*x/kph*M(2*j-2,2*j-4);
    end
    A(j,j) = 1 + 2*(j-1)*(2*j-1)*drot/kph + (drot*(2*j-2)*(2*j-1)+kisc+kph)/(1+drot*(2*j-2)*(2*j-1)+kisc)*x/kph*M(2*j-2,2*j-2);
    if j<nn
        A(j,j+1) = (drot*(2*j-2)*(2*j-1)+kisc+kph)/(1+drot*(2*j-2)*(2*j-1)+kisc)*x/kph*M(2*j-2,2*j);
    end
end
y = A\[1;zeros(nn-1,1)];
b = x*(M(0,0)*y(1)+M(0,2)*y(2))/(1+kisc);
y = 1-y(1)-b;

function y = M(j,k)

if k==j-2
    y = j*(j-1)/(2*j-3)/(2*j-1);
elseif k==j
    y = (2*j*(j+1)-1)/(2*j-1)/(2*j+3);
elseif k==j+2
    y = (j+1)*(j+2)/(2*j+3)/(2*j+5);
end