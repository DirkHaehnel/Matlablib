function y=sigdisp(number,digits)

close all
a=fix(log10(number));
if a>0
if digits>a
    y=number;
else
    y=fix(number/10^(a+1-digits));
end
else
    y=fix(number/10^(a-digits));
end


    
