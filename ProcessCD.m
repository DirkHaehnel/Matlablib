close all
set(0,'defaultaxesfontname','times','defaultaxesfontsize',16);

if 0
    set(0,'defaultaxesfontname','times','defaultaxesfontsize',16);
    pn = 'c:\Temp\';
    dirlist = dir([pn '*.txt']);
    if j>1
        j0=j;
    else
        j0=1;
    end
    for j=j0:length(dirlist)
        line = dirlist(j).name;
        if isempty(findstr(line,'Bl')) & isempty(findstr(line,'Readme'))
            x = load([pn line]);
            if isempty(str2num(line(end-4)))
                bl = load([pn line(1:end-4) 'Bl' line(end-3:end)]);
            else
                bl = load([pn line(1:end-5) 'Bl' line(end-3:end)]);
            end
            s = mean(bl(:,2))-mean(x(:,2));
            if s>0
                [boltzt, boltzx, boltzc] = BoltzmannX(x(:,1), s);   
                line = ['c:/joerg/doc/pampaloni/' line(1:end-4)];
                eval(['print -djpeg -r120 ' line 'Boltzmann']);
                [tx, hx, htx, delay, gaussx, gausstx, kx, cx] = FrancescoX(x(:,1), s);
                eval(['print -djpeg -r120 ' line 'Results']);
                eval(['save ' line 'Results tx hx htx delay gaussx gausstx kx cx s boltzt boltzx boltzc']);
            else
                eval(['save ' line 'NoResult s x bl']);
            end
        end
    end
end

if 1
    namelist = {'Br030801','BrB030801','Br090801','BrC060801','BrB090801','Br060801','BrB060801','BrD060801'};
    conc = [0 10 15 20 25 30 40 50];
    for k=1:length(conc)
        dirlist = dir([namelist{k} '*Results.mat']);
        tmp = [];
        for j=1:length(dirlist)
            line = dirlist(j).name;
            load(line);
            mx = boltzt*boltzx';            
            tmp = [tmp; [cx(2) kx s sum((boltzx-10.^(-polyval(boltzc,boltzt))).^2) sqrt((boltzt-mx).^2*boltzx') boltzc]];
        end
        final{k} = tmp;
    end
end


% com=[]; for j=1:6 com = [com ',t{' num2str(j) '},final{' num2str(j) '}(:,2)']; end; com(1)=[];
% eval(['plot(' com ')']);