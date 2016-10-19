function ExportFigure(name, ftyp)

h = gcf;
P  = get(h,'Position');
set(h,'Position',[P(1) P(2) 900 675]);

h_childs = get(h, 'Children');
if (size(h_childs) == 1)
    TI = get(gca,'TightInset');
    l = TI(1)+0.015;
    r = .97-(TI(1)+TI(3));
    b = TI(2)+0.015;
    if isempty(get(get(gca,'Title'),'string'))
        t = .97-(TI(4)+TI(2));
    else
        t = .935-(TI(4)+TI(2));
    end;
    set(gca,'Position', [l b r t])
end

if nargin<2
    ftyp = 'eps';
end;

if strcmp(ftyp,'eps')
    name = [name '.eps'];
    print('-depsc2', '-r300', name);
end;

if strcmp(ftyp,'png')
    name = [name '.png'];
    print('-dpng', '-r300', '-opengl', name);
end;

if strcmp(ftyp,'tif')
    name = [name '.tif'];
    print('-dtiff', '-r300', name);
end;
