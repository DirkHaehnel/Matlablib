%slowly imlpementation only for bigger structures and out memory error,,
function cum = FASLcumulant0(im,n,flag)

if nargin<3 || isempty(flag)
    im = im - repmat(mean(im,3),[1 1 size(im,3)]);
end

dimensiones = size(im);
switch n
    case 1
        cum = 0;
    case 2
        temporal = zeros(dimensiones(1),dimensiones(2),(dimensiones(3)-1));
%        for j=1:(dimensiones(3)-1) 
%            temporal(:,:,j) = im(:,:,j).*im(:,:,j+1);
%        end
        for f=1:dimensiones(1) 
            temporal(f,:,:) = im(f,:,1:end-1).*im(f,:,2:end);
        end
        cum = mean(temporal,3);
        clear temporal;
    case 3
        temporal = zeros(dimensiones(1),dimensiones(2),(dimensiones(3)-2));
        for f=1:dimensiones(1) 
            temporal(f,:,:) = im(f,:,1:end-2).*im(f,:,2:end-1).*im(f,:,3:end);
        end
        cum = mean(temporal,3);
        clear temporal;
    case 4
%         temporal = zeros(dimensiones(1),dimensiones(2),(dimensiones(3)-3));
%         for f=1:dimensiones(1) 
%             temporal(f,:,:) = im(f,:,1:end-2).*im(f,:,2:end-1).*im(f,:,3:end);
%             
%             im(:,:,1:end-3).*im(:,:,4:end)
%             im(:,:,2:end-2).*im(:,:,3:end-1)
%             
%             im(:,:,1:end-3).*im(:,:,3:end-1)
%             im(:,:,2:end-2).*im(:,:,4:end)
%             im(:,:,1:end-3).*im(:,:,2:end-2)
%             im(:,:,3:end-1).*im(:,:,4:end)
%             im(:,:,1:end-3).*im(:,:,2:end-2).*im(:,:,3:end-1).*im(:,:,4:end)
%         end
%         cum = mean(temporal,3);
%         clear temporal;
%         
        %cum = -mean(,3).*mean(,3)-mean(,3).*mean(,3)-mean(,3).*mean(,3)+mean(,3);
        cum = -mean(im(:,:,1:end-3).*im(:,:,4:end),3).*mean(im(:,:,2:end-2).*im(:,:,3:end-1),3)-mean(im(:,:,1:end-3).*im(:,:,3:end-1),3).*mean(im(:,:,2:end-2).*im(:,:,4:end),3)-mean(im(:,:,1:end-3).*im(:,:,2:end-2),3).*mean(im(:,:,3:end-1).*im(:,:,4:end),3)+mean(im(:,:,1:end-3).*im(:,:,2:end-2).*im(:,:,3:end-1).*im(:,:,4:end),3);
    case 5
        cum = -mean(im(:,:,4:end-1).*im(:,:,5:end),3).*mean(im(:,:,1:end-4).*im(:,:,2:end-3).*im(:,:,3:end-2),3)-mean(im(:,:,3:end-2).*im(:,:,5:end),3).*mean(im(:,:,1:end-4).*im(:,:,2:end-3).*im(:,:,4:end-1),3)-mean(im(:,:,3:end-2).*im(:,:,4:end-1),3).*mean(im(:,:,1:end-4).*im(:,:,2:end-3).*im(:,:,5:end),3)-mean(im(:,:,2:end-3).*im(:,:,5:end),3).*mean(im(:,:,1:end-4).*im(:,:,3:end-2).*im(:,:,4:end-1),3)-mean(im(:,:,2:end-3).*im(:,:,4:end-1),3).*mean(im(:,:,1:end-4).*im(:,:,3:end-2).*im(:,:,5:end),3)-mean(im(:,:,2:end-3).*im(:,:,3:end-2),3).*mean(im(:,:,1:end-4).*im(:,:,4:end-1).*im(:,:,5:end),3)-mean(im(:,:,1:end-4).*im(:,:,5:end),3).*mean(im(:,:,2:end-3).*im(:,:,3:end-2).*im(:,:,4:end-1),3)-mean(im(:,:,1:end-4).*im(:,:,4:end-1),3).*mean(im(:,:,2:end-3).*im(:,:,3:end-2).*im(:,:,5:end),3)-mean(im(:,:,1:end-4).*im(:,:,3:end-2),3).*mean(im(:,:,2:end-3).*im(:,:,4:end-1).*im(:,:,5:end),3)-mean(im(:,:,1:end-4).*im(:,:,2:end-3),3).*mean(im(:,:,3:end-2).*im(:,:,4:end-1).*im(:,:,5:end),3)+mean(im(:,:,1:end-4).*im(:,:,2:end-3).*im(:,:,3:end-2).*im(:,:,4:end-1).*im(:,:,5:end),3);
    case 6
        cum = 2.*mean(im(:,:,1:end-5).*im(:,:,6:end),3).*mean(im(:,:,2:end-4).*im(:,:,5:end-1),3).*mean(im(:,:,3:end-3).*im(:,:,4:end-2),3)+2.*mean(im(:,:,1:end-5).*im(:,:,5:end-1),3).*mean(im(:,:,2:end-4).*im(:,:,6:end),3).*mean(im(:,:,3:end-3).*im(:,:,4:end-2),3)+2.*mean(im(:,:,1:end-5).*im(:,:,6:end),3).*mean(im(:,:,2:end-4).*im(:,:,4:end-2),3).*mean(im(:,:,3:end-3).*im(:,:,5:end-1),3)+2.*mean(im(:,:,1:end-5).*im(:,:,4:end-2),3).*mean(im(:,:,2:end-4).*im(:,:,6:end),3).*mean(im(:,:,3:end-3).*im(:,:,5:end-1),3)+2.*mean(im(:,:,1:end-5).*im(:,:,5:end-1),3).*mean(im(:,:,2:end-4).*im(:,:,4:end-2),3).*mean(im(:,:,3:end-3).*im(:,:,6:end),3)+2.*mean(im(:,:,1:end-5).*im(:,:,4:end-2),3).*mean(im(:,:,2:end-4).*im(:,:,5:end-1),3).*mean(im(:,:,3:end-3).*im(:,:,6:end),3)+2.*mean(im(:,:,1:end-5).*im(:,:,6:end),3).*mean(im(:,:,2:end-4).*im(:,:,3:end-3),3).*mean(im(:,:,4:end-2).*im(:,:,5:end-1),3)+2.*mean(im(:,:,1:end-5).*im(:,:,3:end-3),3).*mean(im(:,:,2:end-4).*im(:,:,6:end),3).*mean(im(:,:,4:end-2).*im(:,:,5:end-1),3)+2.*mean(im(:,:,1:end-5).*im(:,:,2:end-4),3).*mean(im(:,:,3:end-3).*im(:,:,6:end),3).*mean(im(:,:,4:end-2).*im(:,:,5:end-1),3)+2.*mean(im(:,:,1:end-5).*im(:,:,5:end-1),3).*mean(im(:,:,2:end-4).*im(:,:,3:end-3),3).*mean(im(:,:,4:end-2).*im(:,:,6:end),3)+2.*mean(im(:,:,1:end-5).*im(:,:,3:end-3),3).*mean(im(:,:,2:end-4).*im(:,:,5:end-1),3).*mean(im(:,:,4:end-2).*im(:,:,6:end),3)+2.*mean(im(:,:,1:end-5).*im(:,:,2:end-4),3).*mean(im(:,:,3:end-3).*im(:,:,5:end-1),3).*mean(im(:,:,4:end-2).*im(:,:,6:end),3)+2.*mean(im(:,:,1:end-5).*im(:,:,4:end-2),3).*mean(im(:,:,2:end-4).*im(:,:,3:end-3),3).*mean(im(:,:,5:end-1).*im(:,:,6:end),3)+2.*mean(im(:,:,1:end-5).*im(:,:,3:end-3),3).*mean(im(:,:,2:end-4).*im(:,:,4:end-2),3).*mean(im(:,:,5:end-1).*im(:,:,6:end),3)+2.*mean(im(:,:,1:end-5).*im(:,:,2:end-4),3).*mean(im(:,:,3:end-3).*im(:,:,4:end-2),3).*mean(im(:,:,5:end-1).*im(:,:,6:end),3)-mean(im(:,:,1:end-5).*im(:,:,5:end-1).*im(:,:,6:end),3).*mean(im(:,:,2:end-4).*im(:,:,3:end-3).*im(:,:,4:end-2),3)-mean(im(:,:,1:end-5).*im(:,:,4:end-2).*im(:,:,6:end),3).*mean(im(:,:,2:end-4).*im(:,:,3:end-3).*im(:,:,5:end-1),3)-mean(im(:,:,1:end-5).*im(:,:,4:end-2).*im(:,:,5:end-1),3).*mean(im(:,:,2:end-4).*im(:,:,3:end-3).*im(:,:,6:end),3)-mean(im(:,:,1:end-5).*im(:,:,3:end-3).*im(:,:,6:end),3).*mean(im(:,:,2:end-4).*im(:,:,4:end-2).*im(:,:,5:end-1),3)-mean(im(:,:,1:end-5).*im(:,:,3:end-3).*im(:,:,5:end-1),3).*mean(im(:,:,2:end-4).*im(:,:,4:end-2).*im(:,:,6:end),3)-mean(im(:,:,1:end-5).*im(:,:,3:end-3).*im(:,:,4:end-2),3).*mean(im(:,:,2:end-4).*im(:,:,5:end-1).*im(:,:,6:end),3)-mean(im(:,:,1:end-5).*im(:,:,2:end-4).*im(:,:,6:end),3).*mean(im(:,:,3:end-3).*im(:,:,4:end-2).*im(:,:,5:end-1),3)-mean(im(:,:,1:end-5).*im(:,:,2:end-4).*im(:,:,5:end-1),3).*mean(im(:,:,3:end-3).*im(:,:,4:end-2).*im(:,:,6:end),3)-mean(im(:,:,1:end-5).*im(:,:,2:end-4).*im(:,:,4:end-2),3).*mean(im(:,:,3:end-3).*im(:,:,5:end-1).*im(:,:,6:end),3)-mean(im(:,:,1:end-5).*im(:,:,2:end-4).*im(:,:,3:end-3),3).*mean(im(:,:,4:end-2).*im(:,:,5:end-1).*im(:,:,6:end),3)-mean(im(:,:,5:end-1).*im(:,:,6:end),3).*mean(im(:,:,1:end-5).*im(:,:,2:end-4).*im(:,:,3:end-3).*im(:,:,4:end-2),3)-mean(im(:,:,4:end-2).*im(:,:,6:end),3).*mean(im(:,:,1:end-5).*im(:,:,2:end-4).*im(:,:,3:end-3).*im(:,:,5:end-1),3)-mean(im(:,:,4:end-2).*im(:,:,5:end-1),3).*mean(im(:,:,1:end-5).*im(:,:,2:end-4).*im(:,:,3:end-3).*im(:,:,6:end),3)-mean(im(:,:,3:end-3).*im(:,:,6:end),3).*mean(im(:,:,1:end-5).*im(:,:,2:end-4).*im(:,:,4:end-2).*im(:,:,5:end-1),3)-mean(im(:,:,3:end-3).*im(:,:,5:end-1),3).*mean(im(:,:,1:end-5).*im(:,:,2:end-4).*im(:,:,4:end-2).*im(:,:,6:end),3)-mean(im(:,:,3:end-3).*im(:,:,4:end-2),3).*mean(im(:,:,1:end-5).*im(:,:,2:end-4).*im(:,:,5:end-1).*im(:,:,6:end),3)-mean(im(:,:,2:end-4).*im(:,:,6:end),3).*mean(im(:,:,1:end-5).*im(:,:,3:end-3).*im(:,:,4:end-2).*im(:,:,5:end-1),3)-mean(im(:,:,2:end-4).*im(:,:,5:end-1),3).*mean(im(:,:,1:end-5).*im(:,:,3:end-3).*im(:,:,4:end-2).*im(:,:,6:end),3)-mean(im(:,:,2:end-4).*im(:,:,4:end-2),3).*mean(im(:,:,1:end-5).*im(:,:,3:end-3).*im(:,:,5:end-1).*im(:,:,6:end),3)-mean(im(:,:,2:end-4).*im(:,:,3:end-3),3).*mean(im(:,:,1:end-5).*im(:,:,4:end-2).*im(:,:,5:end-1).*im(:,:,6:end),3)-mean(im(:,:,1:end-5).*im(:,:,6:end),3).*mean(im(:,:,2:end-4).*im(:,:,3:end-3).*im(:,:,4:end-2).*im(:,:,5:end-1),3)-mean(im(:,:,1:end-5).*im(:,:,5:end-1),3).*mean(im(:,:,2:end-4).*im(:,:,3:end-3).*im(:,:,4:end-2).*im(:,:,6:end),3)-mean(im(:,:,1:end-5).*im(:,:,4:end-2),3).*mean(im(:,:,2:end-4).*im(:,:,3:end-3).*im(:,:,5:end-1).*im(:,:,6:end),3)-mean(im(:,:,1:end-5).*im(:,:,3:end-3),3).*mean(im(:,:,2:end-4).*im(:,:,4:end-2).*im(:,:,5:end-1).*im(:,:,6:end),3)-mean(im(:,:,1:end-5).*im(:,:,2:end-4),3).*mean(im(:,:,3:end-3).*im(:,:,4:end-2).*im(:,:,5:end-1).*im(:,:,6:end),3)+mean(im(:,:,1:end-5).*im(:,:,2:end-4).*im(:,:,3:end-3).*im(:,:,4:end-2).*im(:,:,5:end-1).*im(:,:,6:end),3);
end


return

z(1)=1; for j=1:1e4 if rand>0.8 z(j+1)=1-z(j); else z(j+1)=z(j); end; end
z = z-mean(z);
zz = permute(z,[3,1,2]);
for j=1:6
    tst(j) = squeeze(cumulant0(zz,j));
end

