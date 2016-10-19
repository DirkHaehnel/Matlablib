function ballplot(x,y,flag,varargin)

if nargin==1
    y = x;
    x = 1:length(x);
end

plot(x,y);
if nargin>3 set(gca,varargin{:}); end
ax = get(gca,'dataaspectratio');
ax = [1 1.7].*ax([1 2]);
fac = 0.05*min(diff(x));
[a,b,c] = ellipsoid(0,0,0,fac,fac,fac*ax(2)/ax(1));
if nargin>2 && ~isstr(flag) && length(flag)==length(x)
    if size(y,1)==1 || size(y,2)==1
        y = y(:);
        surf(x(1)+a,b,y(1)+c,flag(1)*ones(size(a))); 
        hold on; 
        for j=2:length(x) 
            surf(x(j)+a,b,y(j)+c,flag(j)*ones(size(a))); 
        end; 
    else
        x = x(:);
        surf(x(1)+a,b,y(1,1)+c); 
        hold on; 
        for k=1:size(y,2)
            surf(x(1)+a,b,y(1,k)+c,flag(1)*ones(size(a))); 
        end
        for j=2:length(x) 
            for k=1:size(y,2)
                surf(x(j)+a,b,y(j,k)+c,flag(j)*ones(size(a))); 
            end
        end
    end
else
    if size(y,1)==1 || size(y,2)==1
        y = y(:);
        surf(x(1)+a,b,y(1)+c); 
        hold on; 
        for j=2:length(x) 
            surf(x(j)+a,b,y(j)+c); 
        end; 
    else
        x = x(:);
        surf(x(1)+a,b,y(1,1)+c); 
        hold on; 
        for k=1:size(y,2)
            surf(x(1)+a,b,y(1,k)+c); 
        end
        for j=2:length(x) 
            for k=1:size(y,2)
                surf(x(j)+a,b,y(j,k)+c); 
            end
        end
    end
end
hold off
set(gca,'dataaspectratio',ax([1 1 2])); 
shading flat
camlight
set(gca,'ytick',[]);
view([-14, 20]);
set(gcf,'color','w');
set(gca,'color',[1,1,0.9]);
if nargin>3
    for j=1:length(varargin)
        if strcmp(varargin{j}(1),'y') varargin{j}(1)='z'; end
    end
    set(gca,varargin{:});
end
if nargin>2 && isstr(flag) && strcmp(flag,'-')
    hold on 
    plot3(x*ones(1,size(y,2)),x*zeros(1,size(y,2)),y);
    hold off
end
box on
set(gcf,'InvertHardcopy','off');
