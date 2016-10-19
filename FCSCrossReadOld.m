function [y, t] = FCSCrossReadOld(name,t,tcspcfilter)

load(name);

if nargin>1 && ~isempty(t)
    ind = res.autotime>=min(t) & res.autotime<=max(t);
else
    ind = res.autotime>0;
end

if nargin>2 && ~isempty(tcspcfilter)
    y(:,1,:) = res.auto(ind,1,1,:)+res.auto(ind,2,2,:)+res.auto(ind,1,2,:)+res.auto(ind,2,1,:);
    y(:,2,:) = res.auto(ind,3,3,:)+res.auto(ind,4,4,:)+res.auto(ind,3,4,:)+res.auto(ind,4,3,:);
    y(:,3,:) = res.auto(ind,1,4,:)+res.auto(ind,2,3,:)+res.auto(ind,1,3,:)+res.auto(ind,2,4,:);
    y(:,4,:) = res.auto(ind,4,1,:)+res.auto(ind,3,2,:)+res.auto(ind,3,1,:)+res.auto(ind,4,2,:);
else
    y(:,1,:) = sqrt(res.auto(ind,1,2,:).*res.auto(ind,2,1,:));
    y(:,2,:) = sqrt(res.auto(ind,3,4,:).*res.auto(ind,4,3,:));
    y(:,3,:) = sqrt(res.auto(ind,1,4,:).*res.auto(ind,2,3,:));
    y(:,4,:) = sqrt(res.auto(ind,4,1,:).*res.auto(ind,3,2,:));
end

t = res.autotime(ind);
    