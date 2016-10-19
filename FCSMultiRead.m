function [br_br, bb_bb, rr_rr, rb_rb, bb_rr, rr_bb, t] = FCSMultiRead(name,t,filt)

% blue is 0 and red is 1

if ischar(name)
    load(name);
else
    res = name;
end

if nargin>1 && ~isempty(t)
    ind = res.autotime>=min(t) & res.autotime<=max(t);
else
    ind = res.autotime>0;
end

br_br(:,1,:) = res.auto(ind,1,5,:) + res.auto(ind,5,1,:);
br_br(:,2,:) = res.auto(ind,2,6,:) + res.auto(ind,6,2,:);
br_br(:,3,:) = res.auto(ind,1,6,:) + res.auto(ind,5,2,:);
br_br(:,4,:) = res.auto(ind,2,5,:) + res.auto(ind,6,1,:);

bb_bb(:,1,:) = res.auto(ind,9,13,:) + res.auto(ind,13,9,:);
bb_bb(:,2,:) = res.auto(ind,10,14,:) + res.auto(ind,14,10,:);
bb_bb(:,3,:) = res.auto(ind,9,14,:) + res.auto(ind,13,10,:);
bb_bb(:,4,:) = res.auto(ind,10,13,:) + res.auto(ind,14,9,:);

rr_rr(:,1,:) = res.auto(ind,3,7,:) + res.auto(ind,7,3,:);
rr_rr(:,2,:) = res.auto(ind,4,8,:) + res.auto(ind,8,4,:);
rr_rr(:,3,:) = res.auto(ind,3,8,:) + res.auto(ind,7,4,:);
rr_rr(:,4,:) = res.auto(ind,4,7,:) + res.auto(ind,8,3,:);

rb_rb(:,1,:) = res.auto(ind,11,15,:) + res.auto(ind,15,11,:);
rb_rb(:,2,:) = res.auto(ind,12,16,:) + res.auto(ind,16,12,:);
rb_rb(:,3,:) = res.auto(ind,11,16,:) + res.auto(ind,15,12,:);
rb_rb(:,4,:) = res.auto(ind,12,15,:) + res.auto(ind,16,11,:);

rr_bb(:,1,:) = res.auto(ind,3,9,:);
rr_bb(:,2,:) = res.auto(ind,4,10,:);

bb_rr(:,1,:) = res.auto(ind,9,3,:);
bb_rr(:,2,:) = res.auto(ind,10,4,:);


if nargin>2 && ~isempty(filt)
    p = Simplex('ExpFun',3,[],[],[],[],1:size(res.rate,1),sum(res.rate,2));
    [bla, bla, bla, z] = ExpFun(p,1:size(res.rate,1),sum(res.rate,2));
    br_br = br_br(:,:,sum(res.rate,2)<z);
    bb_bb = bb_bb(:,:,sum(res.rate,2)<z);
    rr_rr = rr_rr(:,:,sum(res.rate,2)<z);
    rb_rb = rb_rb(:,:,sum(res.rate,2)<z);
    rr_bb = rr_bb(:,:,sum(res.rate,2)<z);
    bb_rr = bb_rr(:,:,sum(res.rate,2)<z);
end

t = res.autotime(ind);
