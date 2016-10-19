function [tag, tim, head, tch, bin] = ScanRouterRead(name)

% The function [tag, tim, head] = ScanRouterRead(name)
% reads image data from file 'name' for a measurement with two channels in time space
% para(1) = bin time in 100 nsec
% para(2:3) = image size
% tag* = intensity image
% tim* = average lifetime image
% head = measurement parameters

if nargin==0
    [filename, pathname] = uigetfile('*.t3r', 'Interactive mode data:', 0, 0);
    name = [pathname filename];
end

[tmp, tmp, tmp, head] = tttrRead(name,[0 0]);

if head.ScanIdent==1
    [tag, tim, tch, bin] = PIE710Read(name, head);
elseif head.ScanIdent==2
    [tag, tim, tch, bin] = SCX200Read(name);
end

