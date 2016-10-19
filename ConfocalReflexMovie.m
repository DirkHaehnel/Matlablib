rhofield = [0 1.5];
NA = 1.2;
mag = 60;
n0 = 1.333;
nv = 1.52;
dv = 170; % [150:10:190];
lamex = 0.635;
over = [NA*3e3 0 3.5e3];
%atf = [];
atf = 0;
close
for k=1:length(dv)
    %     zimv = -10:0.25:10;
    %     clear tst
    %     for j=1:length(zimv)
    %         [fx0, fx2, by0, by2, rho, phi, intx, inty] = ConfocalReflex([0 0], NA, mag, n0, nv, dv(k), lamex, over, zimv(j), atf);
    %         tst(j) = intx(1);
    %     end
    %
    %     zimv = zimv(tst==max(tst)) + [-0.95:0.1:1];
    zimv = [-0.95:0.1:1] - atf*(1-n0/nv);
    clear intx inty
    for j=1:length(zimv)
        [fx0, fx2, by0, by2, rho, phi, intx(:,:,j), inty(:,:,j)] = ConfocalReflex(rhofield, NA, mag, n0, nv, dv(k), lamex, over, zimv(j), atf);
        subplot(4,5,j)
        pcolor(cos(phi)*rho,sin(phi)*rho,inty(:,:,j)); shading interp; axis image; colormap hot; axis off; %set(gcf,'visible','off')
        title(['defoc = ' num2str(zimv(j)) ' \mum'],'fontsize',10);
        drawnow
        %         if ~isempty(atf)
        %             eval(['print -dpng -zbuffer CF' mint2str(dv(k),3) 'Thickness' mint2str(over(2),2) 'ast' mint2str(atf(2),2) 'cov' mint2str(j,2) 'a']);
        %         else
        %             eval(['print -dpng -zbuffer CF' mint2str(dv(k),3) 'Thickness' mint2str(over(2),2) 'ast00cov' mint2str(j,2) 'a']);
        %         end
    end

    %     if ~isempty(atf)
    %         eval(['save CF' mint2str(dv(k),3) 'Thickness' mint2str(over(2),2) 'ast' mint2str(atf(2),2) 'cov NA mag n0 nv lamex over zimv rho phi intx inty']);
    %     else
    %         eval(['save CF' mint2str(dv(k),3) 'Thickness' mint2str(over(2),2) 'ast00cov NA mag n0 nv lamex over zimv rho phi intx inty']);
    %     end
end

return

for j=1:length(zimv)
    pcolor(cos(phi)*rho,sin(phi)*rho,inty(:,:,j)); shading interp; axis image; colormap gray; axis off
    % caxis([0 max(inty(:))]);
    title(['defoc = ' num2str(zimv(j)) ' \mum']);
    eval(['print -dpng tmp' mint2str(j,2)]);
end

