function im = CombineImages(imraw,n,m,flag,labelx,labely,fsize)

[a,b,c] = size(imraw);
for j=1:n
    for k=1:m
        if (j-1)*m+k<=c
            if nargin>3 && ~isempty(flag) && strcmp(flag,'scale')
                im((j-1)*a+1:j*a,(k-1)*b+1:k*b) = imraw(:,:,(j-1)*m+k)/max(max(imraw(:,:,(j-1)*m+k)));
            else
                im((j-1)*a+1:j*a,(k-1)*b+1:k*b) = imraw(:,:,(j-1)*m+k);
            end
        end
    end
end

mim(im);
hold on
for j=1:n-1
    plot(0.5:size(im,2)+0.5,(j*a+0.5)*ones(1,size(im,2)+1),'y','linewidth',1);
end
for k=1:m-1
    plot((k*b+0.5)*ones(1,size(im,1)+1),0.5:size(im,1)+0.5,'y','linewidth',1);
end
if nargin>4 && ~isempty(labelx) && length(labelx)==m
    for j=1:m
        if nargin>6 && ~isempty(fsize)
            text((2*j-1)*b/2,-size(im,1)/25,labelx{j},'HorizontalAlignment','center','fontsize',fsize)
        else
            text((2*j-1)*b/2,-size(im,1)/25,labelx{j},'HorizontalAlignment','center')
        end
    end
end
if nargin>5 && ~isempty(labely) && length(labely)==n
    for j=1:n
        tmp(j) = length(labely{j});
    end
    tmp = max(tmp);
    for j=1:n
        %text(-5*tmp,(2*j-1)*a/2,labely{j},'HorizontalAlignment','right')
        %text(-size(im,2)/25,(2*j-1)*a/2,labely{j},'HorizontalAlignment','right')
        if nargin>6 && ~isempty(fsize)
            text(-1,(2*j-1)*a/2,labely{j},'HorizontalAlignment','right','fontsize',fsize)
        else
            text(-1,(2*j-1)*a/2,labely{j},'HorizontalAlignment','right')
        end
    end
end
hold off
if nargin>5
    set(gca,'Position',[0.1 0.1 0.9 0.9])
end

if nargout==0
    clear im
end

