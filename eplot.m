function eplot(a,b,varargin)

x{1} = a;
y{1} = b;

cnt = 1;
while cnt<=length(varargin) && isnumeric(varargin{cnt})
    x{(cnt+3)/2} = varargin{cnt};
    y{(cnt+3)/2} = varargin{cnt+1};
    cnt = cnt+2;
end

col ={'r','b','g','m','c','y'};

if iscell(y{1})
    tmp = [];
    for j=1:length(y{1})
        tmp(j) = mean(y{1}{j});
    end
else
    tmp = mean(y{1});
end
plot(x{1},tmp,'o-');
hold on
for j=2:length(x)
    if iscell(y{j})
        tmp = [];
        for k=1:length(y{j})
            tmp(k) = mean(y{j}{k});
        end
    else
        tmp = mean(y{j});
    end
    plot(x{j},tmp,['o-' col{mod(j-1,6)+1}]);
end
ax = axis;
del = diff(ax(1:2))/100;
for k=1:length(x)
    for j=1:length(x{k})
        if iscell(y{k})
            plot(x{k}(j)*[1 1], mean(y{k}{j})+[-1 1]*std(y{k}{j}),col{mod(k-1,6)+1})
            plot(x{k}(j) + [-del del], mean(y{k}{j})+[-1 -1]*std(y{k}{j}),col{mod(k-1,6)+1})
            plot(x{k}(j) + [-del del], mean(y{k}{j})+[1 1]*std(y{k}{j}),col{mod(k-1,6)+1})
        else
            plot(x{k}(j)*[1 1], mean(y{k}(:,j))+[-1 1]*std(y{k}(:,j)),col{mod(k-1,6)+1})
            plot(x{k}(j) + [-del del], mean(y{k}(:,j))+[-1 -1]*std(y{k}(:,j)),col{mod(k-1,6)+1})
            plot(x{k}(j) + [-del del], mean(y{k}(:,j))+[1 1]*std(y{k}(:,j)),col{mod(k-1,6)+1})
        end
    end
end
hold off

