function pd = RotoFit(t, y, pd, name)

% Function RotoFit for fitting Rotational Diffusion data
% (c) Joerg Enderlein, http://www.joerg-enderlein.de (2010)

close all

if (nargin<3) || isempty(pd)
    pd = [100 5e3];
end

expflag = length(pd)-1;
%pd = Simplex('RotoDiffCWFit',pd,[0 pd(2)],[inf pd(2)],[],[], t, y, [1 1]);
pd = Simplex('RotoDiffCWFit',pd,[0 0],[inf inf],[],[], t, y, [1 1]);

if expflag(1)>0
    triplet = pd(end-expflag+1:end);
else
    triplet = [];
end
for j=1:length(pd)-sum(expflag)
    tmp = Rad2RotoDiff(pd(j), 21);
    s{j} = ['\itT\rm_' mint2str(j) ' = ' num2str(tmp,'%.1f') ' ns'];
end
if ~isempty(triplet)
    s{end+1} = '';
    if length(triplet)==1
        s{end+1} = ['\tau = ' num2str(round(triplet(1)/10)/100,'%.1f') ' \mus'];
    else
        tmp = ['\tau = (' num2str(round(triplet(1)/10)/100,'%.1f')];
        for k=2:length(triplet)
            tmp = [tmp ',' num2str(round(triplet(k)/10)/100,'%.1f')];
        end
        tmp = [tmp  ') \mus'];
        s{end+1} = tmp;
    end
end
subplot(4,1,1:3)
text(0.6,0.95,s,'VerticalAlignment','Top','units','normal')
if ~(nargin<4 || isempty(name))
    ExportFigure(name, 'png');
end
    
end


