close all

if 0
    load c:\Joerg\Doc\Microscopy\Hell\Ales\Data2007-11-09\refA4.mat
    load c:\Joerg\Doc\Microscopy\Hell\Ales\Data2007-11-09\fluoA4.mat

    t = (1:2^7)'*0.15;

    plot(t,(sum(sum(ref,3),2)-min(sum(sum(ref,3),2)))/...
        (max(sum(sum(ref,3),2))-min(sum(sum(ref,3),2))),...
        t,(sum(sum(fluo,3),2)-min(sum(sum(fluo,3),2)))/...
        (max(sum(sum(fluo,3),2))-min(sum(sum(fluo,3),2))))

    plot(t,(sum(sum(ref,3),2)-sum(sum(ref(1,:,:),3),2))/...
        (max(sum(sum(ref,3),2))-sum(sum(ref(1,:,:),3),2)),...
        t,(sum(sum(fluo,3),2)-sum(sum(fluo(1,:,:),3),2))/...
        (max(sum(sum(fluo,3),2))-sum(sum(fluo(1,:,:),3),2)) - ...
        t.*(sum(sum(fluo(end,:,:),3),2)-sum(sum(fluo(1,:,:),3),2))/...
        (max(sum(sum(fluo,3),2))-sum(sum(fluo(1,:,:),3),2))/...
        (t(end)-t(1)))

    % close
    % for j=1:60
    %     for k=1:60
    %         plot(t,(sum(sum(ref,3),2)-sum(sum(ref(1,:,:),3),2))/...
    %             (max(sum(sum(ref,3),2))-sum(sum(ref(1,:,:),3),2)),...
    %             t,(fluo(:,j,k)-fluo(1,j,k))/(max(fluo(:,j,k))-fluo(1,j,k)))
    %         drawnow
    %     end
    % end

    %close; p = Simplex('ExpFun',[5 20],[],[],[],[],1:80,sum(sum(fluo(33:100,:,:),3),2),1);
    
    t1 = 33:40;
    t2 = 41:55;
    t3 = 56:75;
    t4 = 76:100;

    %fluo = fluo - repmat(fluo(100,:,:),[size(fluo,1),1,1]);
    fluo = cat(1,sum(fluo(t1,:,:)),sum(fluo(t2,:,:)),sum(fluo(t3,:,:)),sum(fluo(t4,:,:)));

    t = 1:78;
    rate = 1./[5 50];
    master = [exp(-t'*rate) ones(size(t'))];
    master = master./(ones(size(t'))*sum(master));
    master = [sum(master(t1-32,:)); sum(master(t2-32,:)); sum(master(t3-32,:)); sum(master(t4-32,:))];
    
    clear para;
    for j=1:60
        for k=1:60
            para(:,j,k) = lsqnonneg(master,fluo(:,j,k));
            %para(:,j,k) = Simplex('ExpFun', 2, 10, inf, [], [], t, fluo(:,j,k));
            %         plot(t,fluo(:,j,k),'o',t,master*para(:,j,k));
            %         drawnow
        end
    end
    mim(cat(3,squeeze(sum(fluo)),mConv2(squeeze(para(1,:,:)),disk(1))))
end

if 1
    load c:\Joerg\Doc\Microscopy\Hell\Ales\Data2007-11-09\Scan4.mat
    fluo0 = ScanSmall4;

    %close all; p = Simplex('ExpFun',[5 20],[],[],[],[],1:80,sum(sum(fluo(21:100,:,:),3),2),1);

    t1 = 21:25;
    t2 = 26:35;
    t3 = 36:55;
    t4 = 56:100;

    fluo = fluo0 - repmat(fluo0(100,:,:),[size(fluo0,1),1,1]);
    fluo = cat(1,sum(fluo(t1,:,:)),sum(fluo(t2,:,:)),sum(fluo(t3,:,:)),sum(fluo(t4,:,:)));

    t = 1:80;
    rate = 1./[5 40];
    %master = [exp(-(t-21)'*rate) ones(size(t'))];
    master = exp(-t'*rate);
    master = master./(ones(size(t'))*sum(master));
    master = [sum(master(t1-20,:)); sum(master(t2-20,:)); sum(master(t3-20,:)); sum(master(t4-20,:))];
    
    clear para c
    for j=1:50
        for k=1:50
            %para(:,j,k) = master\fluo(:,j,k);
            para(:,j,k) = lsqnonneg(master,fluo(:,j,k));
            %para(:,j,k) = sort(Simplex('ExpFun', [10 50], [0 0], [], [], [], t, fluo(t,j,k)));
            %[err, c(:,j,k)] = ExpFun(para(:,j,k), t, fluo(t,j,k), 1);
            %         plot(t,fluo(:,j,k),'o',t,master*para(:,j,k));
            %         drawnow
        end
    end
    %mim(squeeze(para(1,:,:)))
    
    mim(squeeze(sum(fluo0)),squeeze(sum(fluo0))); 
    figure; 
    mim(mConv2(squeeze(para(1,:,:)),disk(2)),squeeze(sum(fluo0)))
    
    %mim(squeeze(para(1,:,:)),squeeze(sum(fluo)))
end