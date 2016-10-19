load c:\Joerg\Doc\Fcs\Samantha\SamData.mat

anti1 = mHist(diff(x(4,:)),0:5000:1e6,double(x(2,1:end-1)==0).*double(x(2,2:end)==1));
[anti2 t] = mHist(diff(x(4,:)),0:5000:1e6,double(x(2,2:end)==0).*double(x(2,1:end-1)==1));

plot(4e-3*[-fliplr(t) t],[flipud(anti2); anti1])
xlabel('time (ns)')
ylabel('autocorrelation (a.u.)');