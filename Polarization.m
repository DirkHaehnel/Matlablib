namea = 'c:\Joerg\Doc\Boehmer\290501\3DOrientierung\2Ja167\3AJa167.t3r';
nameb = 'c:\Joerg\Doc\Boehmer\290501\3DOrientierung\2Ja167\3BJa167.t3r';
para = [5000 500 500];

[tag1a, tag2a, tim1a, tim2a, rela, thist1a, thist2a, va, head] = Scan2Read(namea, para);
[tag1b, tag2b, tim1b, tim2b, relb, thist1b, thist2b, vb, head] = Scan2Read(nameb, para);

[ta1a, ta2a, ta1b, ta2b, ti1a, ti2a, ti1b, ti2b, c] = PolPicAlign(tag1a, tag2a, tag1b, tag2b, tim1a, tim2a, tim1b, tim2b);

%ta2a = max(max(ta1a))/max(max(ta2a))*ta2a;
%ta2b = max(max(ta1b))/max(max(ta2b))*ta2b;

%ta2a = ta2a/rela;
%ta2b = ta2b/relb;

phi = atan2(ta1a-ta2a+ta1b-ta2b,ta1a+ta2a-ta1b-ta2b)/2;
para = sqrt((ta1a-ta2a+ta1b-ta2b).^2 + (ta1a+ta2a-ta1b-ta2b).^2);
ortho = abs((ta1a+ta2a+ta1b+ta2b-para)/4);
% for j=1:50 mimArrow(para,phi,(j-1)*pi/25); eval(['print -djpeg -r100 tmp' num2str(floor(j/10)) num2str(rem(j,10))]); end
% c:/imagemagick/convert +map tmp*.jpg tst.gif

tag = ta1a+ta2a+ta1b+ta2b;
tim = (ti1a.*ta1a+ti2a.*ta2a+ti1b.*ta1b+ti2b.*ta2b)./tag;
[mol, len, phot, int, tav, rim] = DataPro2D(tag, tim);

phiv = zeros(size(len));
parav = zeros(size(len));
orthov = zeros(size(len));
for j=1:length(len)
   ind = (mol==j);
   cp = sum(ta1a(ind)-ta2a(ind)+ta1b(ind)-ta2b(ind));
   sp = sum(ta1a(ind)+ta2a(ind)-ta1b(ind)-ta2b(ind));
   parav(j) = sqrt(sp^2 + cp^2);
   phiv(j) = atan2(cp,sp);
   orthov(j) = (sum(ta1a(ind)+ta2a(ind)+ta1b(ind)+ta2b(ind))-parav(j))/4;
end    