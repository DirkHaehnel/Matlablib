% IrisMaster

direx = {'D:\Joerg\Doc\Iris\061103_TMR-DOPE_SPB_Weitfeld'};

if 1
    for kk=1:length(direx)
        names = dir([direx{kk} '*4.tif']);
        for jj=1:length(names)
            name = [direx{kk} names(jj).name];

            [im, cim, resx, resy, im0] = Iris(name);

            imlen = size(im,3);

            cx = [min(im0(:)) max(im0(:))];
            for j=1:imlen
                mim(im0(:,:,j));
                caxis(cx);
                mov(j) = getframe;
            end
            movie2avi(mov,[name '_raw.avi'],'compression','none');

            cx = [min(im(:)) max(im(:))];
            for j=1:imlen
                mim(cim(:,:,j));
                caxis(cx);
                hold on;
                plot(resx{j},resy{j},'ob');
                hold off
                mov(j) = getframe;
            end
            movie2avi(mov,[name '.avi'],'compression','none');

            bin = 0:ceil(size(im,1)^2+size(im,2).^2);
            dist = zeros(length(bin),imlen-1);
            for k=1:imlen-1
                cnt = 0;
                for j=1:imlen-k
                    tmp = (resx{j}'*ones(1,length(resx{j+k}))-ones(length(resx{j}),1)*resx{j+k}).^2 + ...
                        (resy{j}'*ones(1,length(resy{j+k}))-ones(length(resy{j}),1)*resy{j+k}).^2;
                    dist(:,k) = dist(:,k) + mHist(tmp(:),bin);
                    cnt = cnt + prod(size(tmp));
                end
                dist(:,k) = dist(:,k)/cnt;
            end
            c = polyfit(1:imlen-1,bin*dist./sum(dist),1);
            plot(1:imlen-1,bin*dist./sum(dist),1:imlen-1,polyval(c,1:imlen-1));
            xlabel('frame #');
            ylabel('\langle\itr\rm^2\rangle [pix^2]')
            title(['slope = ' mnum2str(c(1),2,1) ' pix^2/T_{frame}'])
            eval(['print -dpng ''' name '''.png' ])

        end
    end
end

res = [];
if 0
    for kk=1:3
        names = dir([direx{kk} '20*']);
        for jj=1:length(names)
            if ~strcmp(names(jj).name(end-2:end),'gif') & ~strcmp(names(jj).name(end-2:end),'png')
                name = [direx{kk} names(jj).name];

                [im, cim, resx, resy, im0, bck] = Iris(name);

                imlen = size(im,3);

                auto = zeros(round(imlen/2),3);

                for j=1:round(imlen/2)
                    for k=1:imlen-j
                        auto(j,1) = auto(j,1) + sum(sum(im0(:,:,k).*im0(:,:,k+j)));
                        auto(j,2) = auto(j,2) + sum(sum(im(:,:,k).*im(:,:,k+j)));
                        auto(j,3) = auto(j,3) + sum(sum(cim(:,:,k).*cim(:,:,k+j)));
                    end;
                end

                auto = auto./((imlen-1:-1:imlen-size(auto,1))'*ones(1,size(auto,2)));

                plot(1:size(auto,1),auto./(ones(size(auto,1),1)*max(auto)))
                xlabel('frame #');
                ylabel('autocorrelation')
%                eval(['print -dpng ''' name '''_auto.png' ])

                res = [res simplex('rigler',[10 0],[0 0],[inf eps],[],[],1:size(auto,1),auto(:,3))];
            end
        end
    end
end

