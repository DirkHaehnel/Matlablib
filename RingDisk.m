function z = ring(r1,r2,l,b,x,y);

z = zeros(l,b);

for a=1:l
    for b=1:b
        t = ((x-a)^2+(y-b)^2)^(1/2);
        t1 = (r2-t);
        t2 = (t-r1);
        if (t1>=0)&(t2>=0)
            z(a,b) = 1;
        end;
    end;
end;

t = sum(sum(z));
if t>0
    z = z./t;
end;