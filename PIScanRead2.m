function [intens, flim, tcspc, head] = PIScanRead2(name)

[y, tcspc, chan, markers, num, overcount, head] = pt3v2read(name, [1 inf]);

intens = zeros(head.ImgHdr.PixY,head.ImgHdr.PixX);

Turns = y(markers==4 & chan==15);

h = waitbar(0,'Reading image data');

if head.ImgHdr.Pattern==0
    for j=1:length(Turns)-1
        dT = (Turns(j+1)-Turns(j)); 
        t1 = Turns(j)+head.ImgHdr.TStartTo*dT;
        t2 = Turns(j)+head.ImgHdr.TStopTo*dT;
        tmp = y(y>=t1 & y<=t2) - t1;
        intens(j,:) = tttr2bin([0; tmp; t2-t1-1], (t2-t1)/(head.ImgHdr.PixX));
        intens(j,[1 end]) = intens(j,[1 end]) - 1;
        tmp = tcspc(y>=t1 & y<=t2);
        cnt = cumsum([1 intens(j,:)]);
        for k=1:head.ImgHdr.PixX
            if cnt(k+1)>cnt(k)
                flim{j+1,head.ImgHdr.PixX-k+1} = tmp(cnt(k):cnt(k+1)-1);
            else
                flim{j+1,k} = 0;
            end
        end
        intens(j+1,:) = intens(j+1,end);
        waitbar(j/length(Turns),h); drawnow
    end
else    
    for j=1:2:length(Turns)-1
        dT = (Turns(j+1)-Turns(j));
        t1 = Turns(j)+head.ImgHdr.TStartTo*dT;
        t2 = Turns(j)+head.ImgHdr.TStopTo*dT;
        tmp = y(y>=t1 & y<=t2) - t1;
        intens(j,:) = tttr2bin([0; tmp; t2-t1-1], (t2-t1)/(head.ImgHdr.PixX));
        intens(j,[1 end]) = intens(j,[1 end]) - 1;
        tmp = tcspc(y>=t1 & y<=t2);
        cnt = cumsum([1 intens(j,:)]);
        for k=1:head.ImgHdr.PixX
            if cnt(k+1)>cnt(k)
                flim{j,k} = tmp(cnt(k):cnt(k+1)-1);
            else
                flim{j,k} = 0;
            end
        end
        t1 = Turns(j)+head.ImgHdr.TStartFro*dT;
        t2 = Turns(j)+head.ImgHdr.TStopFro*dT;
        tmp = y(y>=t1 & y<=t2) - t1;
        intens(j+1,:) = tttr2bin([0; tmp; t2-t1-1], (t2-t1)/(head.ImgHdr.PixX));
        intens(j+1,[1 end]) = intens(j+1,[1 end]) - 1;
        tmp = tcspc(y>=t1 & y<=t2);
        cnt = cumsum([1 intens(j+1,:)]);
        for k=1:head.ImgHdr.PixX
            if cnt(k+1)>cnt(k)
                flim{j+1,head.ImgHdr.PixX-k+1} = tmp(cnt(k):cnt(k+1)-1);
            else
                flim{j+1,k} = 0;
            end
        end
        intens(j+1,:) = intens(j+1,end:-1:1);
        waitbar(j/length(Turns),h); drawnow
    end
end
tcspc = mHist(tcspc,0:2^12-1);
close(h)


