x=[1	0.7	4.4	0.1	0.1	
2	0.8	3.9	0.1	0.2
3	0.6	4.0	0.1	0.2
4	0.7	4.1	0.1	0.1
5	0.7	4.2	0.2	0.2
6	0.6	4.2	0.1	0.1
7	0.6	4.1	0.1	0.2
8	0.7	4.5	0.2	0.2
9	0.7	4.1	0.1	0.1
10	0.8	4.5	0.1	0.1
11	0.8	4.2	0.1	0.1
12	0.8	4.2	0.1	0.2
13	0.7	4.2	0.2	0.2
14	0.9	4.2	0.1	0.1
15	0.7	3.7	0.1	0.1
16	0.8	4.1	0.1	0.1
17	0.6	4.1	0.1	0.2
18	0.8	4.3	0.2	0.2
19	0.7	3.8	0.2	0.2

20	0.5	4.0	0.1	0.2
21	0.8	4.1	0.1	0.1
22	0.7	4.0	0.1	0.1
23	0.7	3.9	0.1	0.1
24	0.8	4.3	0.2	0.2
25	0.6	3.9	0.1	0.2
26	0.7	4.1	0.1	0.2
27	0.7	4.3	0.1	0.1
28	0.7	4.0	0.1	0.2];

[xx,yy]=meshgrid(0.4:0.01:1,3.4:0.01:5);
tst=0*xx;
for k=1:28 
    tst = tst + exp(-(xx-x(k,2)).^2/2/x(k,4)^2 - (yy-x(k,3)).^2/2/x(k,5)^2)/(2*pi*x(k,4)*x(k,5)); 
end; 
tst = tst/28;
pcolor(xx,yy,tst)
shading interp
colormap hot
xlabel('quantum yield'); 
ylabel('free lifetime (ns)')

c = lscov(xx(:),yy(:),tst(:));
hold on 
plot(xx(:),1/c*(xx(:)-sum(sum(tst.*xx))/sum(tst(:)))+sum(sum(tst.*yy))/sum(tst(:)),'c')
plot(xx(1,:),ones(1,size(xx,2))*sum(sum(tst.*yy))/sum(tst(:)),'b')
hold off