% Atto655-Trp data analysis

if 0
    cd D:\Daten\Gregor\2005-09-14_Piehler
    names = dir('*.pt3');
    for j=1:2 [anti(j) head(j)]=pt3pro(names(2*j-1).name,names(2*j).name,'antibunch',[3125 500]); end
end

if 1
   close
   p1=simplex('expfun',1e3,[0],[],[],[],anti(1).autotime,mean(anti(1).auto,2),1); 
   p2=simplex('expfun',1e3,[0],[],[],[],anti(2).autotime,mean(anti(2).auto,2),1);
   
   [err1, c1, zz1, z1] = ExpFun(p1, anti(1).autotime,mean(anti(1).auto,2));
   [err2, c2, zz2, z2] = ExpFun(p2, anti(2).autotime,mean(anti(2).auto,2));

   plot(anti(1).autotime,z1,'r',anti(2).autotime,z2,'b',anti(1).autotime,mean(anti(1).auto,2),'or',anti(2).autotime,mean(anti(2).auto,2),'ob')
   %    plot(anti(1).autotime,z1/max(z1),'r',anti(2).autotime,z2/max(z2),'b',...
   %        anti(1).autotime,mean(anti(1).auto,2)/max(z1),'or',anti(2).autotime,mean(anti(2).auto,2)/max(z2),'ob')


end

