function y = BeadScanRead0(fname)

read_fun = 'SCX200Read';
read_fun = 'PIE710Read';

close all hidden

eval(['[tag, tim, tch, bin] = ' read_fun '(fname);'])

semilogy(bin,tch);
pos = ginput(2);

filter = bin>pos(1,1) & bin<pos(2,1);
eval(['[tag1, tim1] = ' read_fun '(fname,[filter'' filter'']);'])
filter = bin<pos(1,1) | bin>pos(2,1);
eval(['[tag2, tim2] = ' read_fun '(fname,[filter'' filter'']);'])

y(:,:,1) = sum(tim1,3);
y(:,:,2) = sum(tim2,3);