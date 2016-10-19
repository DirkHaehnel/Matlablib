name = 'c:\Joerg\Doc\Fcs\BSA\antibunching_bsa_atto647_without_dic_189C.pt3';

photons = 1e6;
[y,x,num,num,num] = pt3v2read(name, [1 photons]);
[b a] = mHist(x, min(x):max(x));
if num == photons 
    tst=1; 
    cnts = num;
else
    tst=0; 
end
while tst
    [y,x,num,num,num] = pt3v2read(name, [cnts+1 photons]);
    b = b + mHist(x,a);
    cnts = cnts + num;
    if num<photons
        tst = 0;
    end
    plot(a(2:end)*32*1e-3,b(2:end))
    drawnow
end

save antibunching a b

