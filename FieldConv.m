function fout = FieldConv(f1,f2)
    
maxm = (size(f1,3)-1)/2;
fout = 0*f1;

fout(:,:,1,:) = f1(:,:,1,:).*f2(:,:,1,:);
for j=1:maxm
    fout(:,:,1+j,:) = f1(:,:,1,:).*f2(:,:,1+j,:) + f2(:,:,1,:).*f1(:,:,1+j,:);
    fout(:,:,maxm+1+j,:) = f1(:,:,1,:).*f2(:,:,maxm+1+j,:) + f2(:,:,1,:).*f1(:,:,maxm+1+j,:);
end

for j1=1:maxm
    for j2=1:maxm
        tmp = f1(:,:,1+j1,:).*f2(:,:,1+j2,:);
        fout(:,:,1+abs(j1-j2),:) = fout(:,:,1+abs(j1-j2),:) + 0.5*tmp;
        if j1+j2<=maxm
            fout(:,:,1+j1+j2,:) = fout(:,:,1+j1+j2,:) + 0.5*tmp;
        end
        tmp = f1(:,:,maxm+1+j1).*f2(:,:,maxm+1+j2);
        fout(:,:,1+abs(j1-j2),:) = fout(:,:,1+abs(j1-j2),:) - 0.5*tmp;
        if j1+j2<=maxm
            fout(:,:,1+j1+j2,:) = fout(:,:,1+j1+j2,:) + 0.5*tmp;
        end
        tmp = f1(:,:,1+j1,:).*f2(:,:,maxm+1+j2);
        fout(:,:,1+maxm+abs(j1-j2),:) = fout(:,:,maxm+1+abs(j1-j2),:) + 0.5*tmp;
        if j1+j2<=maxm
            fout(:,:,maxm+1+j1+j2,:) = fout(:,:,maxm+1+j1+j2,:) + 0.5*tmp;
        end
        tmp = f1(:,:,maxm+1+j1,:).*f2(:,:,1+j2,:);
        fout(:,:,1+maxm+abs(j1-j2),:) = fout(:,:,maxm+1+abs(j1-j2),:) + 0.5*tmp;
        if j1+j2<=maxm
            fout(:,:,maxm+1+j1+j2,:) = fout(:,:,maxm+1+j1+j2,:) + 0.5*tmp;
        end
    end
end
