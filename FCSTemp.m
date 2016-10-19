close all
%     t = res(1).autotime;
%     for j=1:length(res)
%         x(:,j) = mean(res(j).auto')';
%     end
%     tt = 10:length(t);
pow = [0.2000    0.6670    2.0000    6.6700   20.0000   66.7000]*5;    
jo0 = jo;
jz0 = jz;
ja0 = ja;
js0 = js; 
for jo=jo0:size(modres,2)
    if jo==jo0 jz1 = jz0; else jz1 = 1; end
    for jz=jz1:2; %size(modres,3)
        if (jo==jo0 & jz==jz0) ja1 = ja0; else ja1 = 1; end            
        for ja=ja1:size(modres,4)
            if (jo==jo0 & jz==jz0 & ja==ja0) js1 = js0; else js1 = 1; end                            
            for js=js1:size(modres,5)
                p(:,jo,jz,ja,js) = simplex('alexafit',[5e-7 1e5 1e3],zeros(1,3),[],[],[],pow,t(tt),x(tt,6:-1:1),autotime,squeeze(modres(:,jo,jz,ja,js,1:6)));
%                p(:,jo,jz,ja,js) = simplex('alexafit',p(:,jo,jz,ja,js),zeros(1,3),[],[],[],pow,t(tt),x(tt,6:-1:1),autotime,squeeze(modres(:,jo,jz,ja,js,1:6)));
                [err(jo,jz,ja,js), bla, z(:,jo,jz,ja,js,:)] = alexafit(p(:,jo,jz,ja,js),pow,t(tt),x(tt,6:-1:1),autotime,squeeze(modres(:,jo,jz,ja,js,1:6)));
            end
        end
    end
end

tst = err(1,1,1,1);
ind = [1 1 1 1];
jo0 = 1;
jz0 = 1;
ja0 = 1;
js0 = 1; 
for jo=jo0:size(modres,2)
    if jo==jo0 jz1 = jz0; else jz1 = 1; end
    for jz=jz1:2; %size(modres,3)
        if (jo==jo0 & jz==jz0) ja1 = ja0; else ja1 = 1; end            
        for ja=ja1:size(modres,4)
            if (jo==jo0 & jz==jz0 & ja==ja0) js1 = js0; else js1 = 1; end                            
            for js=js1:size(modres,5)
                if err(jo,jz,ja,js)<tst
                    tst = err(jo,jz,ja,js);
                    ind = [jo,jz,ja,js];
                end
            end
        end
    end
end

    
