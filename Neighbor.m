function [pthx, pthy, len, trck, x, y] = Neighbor(x,y)

% (c) Joerg Enderlein (2006)

if iscell(x)
    xx(:,1) = x{1}(:);
    yy(:,1) = y{1}(:);
    for j=2:length(x)
        mm = size(xx,1);
        len = length(x{j});
        if len>mm
            xx(mm+1:len,:) = NaN;
            yy(mm+1:len,:) = NaN;
            xx(:,j) = x{j}(:);
            yy(:,j) = y{j}(:);
        else
            xx(1:len,j) = x{j}(:);
            yy(1:len,j) = y{j}(:);
            xx(len+1:mm,j) = NaN;
            yy(len+1:mm,j) = NaN;
        end
    end
    x = xx;
    y = yy;
end

vec = 1:size(x,1);
col = ones(size(x,1),1);
ngbr = zeros(size(x,1),size(x,2)-1);
for j=1:size(x,2)-1
    dist = (col*x(:,j)'-x(:,j+1)*col').^2+(col*y(:,j)'-y(:,j+1)*col').^2;
    [~, ind] = sort(dist(:));
    ind1 = col*vec;
    ind2 = vec'*col';
    ind1 = ind1(ind);
    ind2 = ind2(ind);
    ind = 0*col;
    while ~isempty(ind1)
        k = ind1(1);
        ind(k) = ind2(1);
        ind1(1) = []; ind2(1) = [];
        if ~isempty(ind1) && sum(ind2==ind(k))>0
            ind1(ind2==ind(k)) = [];
            ind2(ind2==ind(k)) = [];
        end
        if ~isempty(ind1) && sum(ind1==k)>0
            ind2(ind1==k) = [];
            ind1(ind1==k) = [];
        end
    end
    ngbr(isfinite(x(:,j)) & isfinite(x(ind,j+1)),j) = ind(isfinite(x(:,j)) & isfinite(x(ind,j+1)));
end
trck = ngbr;

if nargout>1
    col = 1;
    cnt = 1;
    while col<=size(ngbr,2)
        if all(ngbr(:,col)==0)
            col = col+1;
        else
            j = 1;
            while ngbr(j,col)==0 
                j = j+1; 
            end
            k = col;
            pthx{cnt}(k-col+1) = x(j,k);            
            pthy{cnt}(k-col+1) = y(j,k);
            while k<=size(ngbr,2) && ngbr(j,k)>0
                tmp = ngbr(j,k); ngbr(j,k) = 0; j = tmp;
                k = k+1;
                pthx{cnt}(k-col+1) = x(j,k);            
                pthy{cnt}(k-col+1) = y(j,k);
            end
            cnt = cnt+1;
        end
    end
    
    for j=1:length(pthx)
        len(j) = length(pthx{j});
    end
    pthx(len==1)=[];
    pthy(len==1)=[];
    
    len = length(pthx);
    for j=1:len
        tmp = [tmp sqrt(diff(pthx{j}).^2+diff(pthy{j}).^2)];
    end
    tmp = 0.8*mean(tmp);
    for j=1:len
        ind = 1:length(pthx{j})-1;
        ind = [ind(sqrt(diff(pthx{j}).^2+diff(pthy{j}).^2)>tmp) length(pthx{j})];
        for k=1:length(ind)-1
            pthx{end} = pthx{j}(ind(k)+1:ind(k+1));
            pthy{end} = pthy{j}(ind(k)+1:ind(k+1));
        end
        pthx{j} = pthx{j}(1:ind(1));
        pthy{j} = pthy{j}(1:ind(1));        
    end
    for j=1:length(pthx)
        len(j) = length(pthx{j});
    end

    pthx(len==1)=[];
    pthy(len==1)=[];
    len(len==1)=[];
    [len,ind]=sort(len);
    pthx=pthx(ind);
    pthy=pthy(ind);
    
    %frb = ['mcrgb'];
    plot(pthx{1},pthy{1});%,frb(1));
    hold on;
    for j=2:length(pthx)
        %plot(pthx{j},pthy{j},frb(mod(j-1,5)+1));
        plot(pthx{j},pthy{j});
    end;
    hold off
    axis image
    colorize
    
end
                    
