%load CavityAbsorption

pcolor([outer_radius_vec max(outer_radius_vec)-min(outer_radius_vec)+outer_radius_vec],rm,[log10(emx) log10(emz)]); 
xlabel('outer radius (nm)');
ylabel('emitter position (nm)');
axis image; 
shading interp; 
caxis([-4 6]); 
set(gca,'linewidth',0.5)
patch([2.5, 5, max(outer_radius_vec), max(outer_radius_vec), 2.5], [0, 0, max(outer_radius_vec)-5, max(outer_radius_vec), 2.5],[255 250 250]/255); 
patch(max(outer_radius_vec)-min(outer_radius_vec)+[2.5, 5, max(outer_radius_vec), max(outer_radius_vec), 2.5], ...
    [0, 0, max(outer_radius_vec)-5, max(outer_radius_vec), 2.5],[255 250 250]/255); 
line(max(outer_radius_vec)*[1 1],[min(rm) max(rm)],'color','b','linewidth',1)
set(gca,'xtick',[5:5:25 max(outer_radius_vec)-min(outer_radius_vec)+(5:5:25)])
set(gca,'xticklabel',[5:5:25 5:5:25])
text(0.5*max(outer_radius_vec)-0.*min(outer_radius_vec),42.5,'lg {\it S}_{\rm ||} /{\it S}_{\rm 0||}','horizontalalignment','center')
text(1.5*max(outer_radius_vec)-min(outer_radius_vec),42.5,'lg {\it S}_{\rm \perp} /{\it S}_{\rm 0\perp}','horizontalalignment','center')
h = colorbar;
%set(get(h,'ylabel'),'string','lg \it{S}','interpreter','tex')
set(h,'linewidth',0.5)
print -djpeg -r300 CavityAbsorptionS


pcolor([outer_radius_vec max(outer_radius_vec)-min(outer_radius_vec)+outer_radius_vec],rm,[log10(intx) log10(intz)]); 
xlabel('outer radius (nm)');
ylabel('emitter position (nm)');
axis image; 
shading interp; 
caxis([-4 6]); 
set(gca,'linewidth',0.5)
patch([2.5, 5, max(outer_radius_vec), max(outer_radius_vec), 2.5], [0, 0, max(outer_radius_vec)-5, max(outer_radius_vec), 2.5],[255 250 250]/255); 
patch(max(outer_radius_vec)-min(outer_radius_vec)+[2.5, 5, max(outer_radius_vec), max(outer_radius_vec), 2.5], ...
    [0, 0, max(outer_radius_vec)-5, max(outer_radius_vec), 2.5],[255 250 250]/255); 
line(max(outer_radius_vec)*[1 1],[min(rm) max(rm)],'color','b','linewidth',1)
set(gca,'xtick',[5:5:25 max(outer_radius_vec)-min(outer_radius_vec)+(5:5:25)])
set(gca,'xticklabel',[5:5:25 5:5:25])
text(0.5*max(outer_radius_vec)-0.*min(outer_radius_vec),42.5,'lg \langle|\bf{E}\rm^{\it{sc}}_{\rm ||}|^2\rangle /\langle|\bf{E}\rm_{\rm ||}|^2\rangle_0','horizontalalignment','center')
text(1.5*max(outer_radius_vec)-min(outer_radius_vec),42.5,'lg \langle|\bf{E}\rm^{\it{sc}}_{\rm \perp}|^2\rangle /\langle|\bf{E}\rm_{\rm \perp}|^2\rangle_0','horizontalalignment','center')
h = colorbar;
%set(get(h,'ylabel'),'string','lg \langle|\bf{E}\rm^{\it{sc}}|^2\rangle')
set(h,'linewidth',0.5)
print -djpeg -r300 CavityAbsorptionEsc


pcolor([outer_radius_vec max(outer_radius_vec)-min(outer_radius_vec)+outer_radius_vec],rm,[log10(flx) log10(flz)]); 
xlabel('outer radius (nm)');
ylabel('emitter position (nm)');
axis image; 
shading interp; 
caxis([-4 6]); 
set(gca,'linewidth',0.5)
patch([2.5, 5, max(outer_radius_vec), max(outer_radius_vec), 2.5], [0, 0, max(outer_radius_vec)-5, max(outer_radius_vec), 2.5],[255 250 250]/255); 
patch(max(outer_radius_vec)-min(outer_radius_vec)+[2.5, 5, max(outer_radius_vec), max(outer_radius_vec), 2.5], ...
    [0, 0, max(outer_radius_vec)-5, max(outer_radius_vec), 2.5],[255 250 250]/255); 
line(max(outer_radius_vec)*[1 1],[min(rm) max(rm)],'color','b','linewidth',1)
set(gca,'xtick',[5:5:25 max(outer_radius_vec)-min(outer_radius_vec)+(5:5:25)])
set(gca,'xticklabel',[5:5:25 5:5:25])
text(0.5*max(outer_radius_vec)-0.*min(outer_radius_vec),42.5,'lg \langle|\bf{E}\rm^{\it{re}}_{\rm ||}|^2\rangle /\langle|\bf{E}\rm_{\rm ||}|^2\rangle_0','horizontalalignment','center')
text(1.5*max(outer_radius_vec)-min(outer_radius_vec),42.5,'lg \langle|\bf{E}\rm^{\it{re}}_{\rm \perp}|^2\rangle /\langle|\bf{E}\rm_{\rm \perp}|^2\rangle_0','horizontalalignment','center')
h = colorbar;
%set(get(h,'ylabel'),'string','lg \langle|\bf{E}\rm^{\it{re}}|^2\rangle')
set(h,'linewidth',0.5)
print -djpeg -r300 CavityAbsorptionEre



% pcolor(outer_radius_vec,rm,log10(emz)); 
% xlabel('outer radius (nm)');
% ylabel('emitter position (nm)');
% axis image; 
% shading interp; 
% caxis([-4 6]); 
% patch([2.5, 5, 30, 30, 2.5], [0, 0, 25, 30, 2.5],[255 250 250]/255); 
% h = colorbar;
% set(get(h,'ylabel'),'string','lg {\itS}_{\rm \perp}','interpreter','tex')
% print -djpeg -r300 CavityAbsorptionSz

% pcolor(outer_radius_vec,rm,log10(intz)); 
% xlabel('outer radius (nm)');
% ylabel('emitter position (nm)');
% axis image; 
% shading interp; 
% caxis([-4 6])
% patch([2.5, 5, 30, 30, 2.5], [0, 0, 25, 30, 2.5],[255 250 250]/255); 
% h = colorbar;
% set(get(h,'ylabel'),'string','lg |\bf{E}\rm^{\it{sc}}_{\rm \perp}|^2','interpreter','tex')
% print -djpeg -r300 CavityAbsorptionEscz

% pcolor(outer_radius_vec,rm,log10(flz)); 
% xlabel('outer radius (nm)');
% ylabel('emitter position (nm)');
% axis image; 
% shading interp; 
% caxis([-4 6]); 
% patch([2.5, 5, 30, 30, 2.5], [0, 0, 25, 30, 2.5],[255 250 250]/255); 
% h = colorbar;
% set(get(h,'ylabel'),'string','lg |\bf{E}\rm^{\it{re}}_{\rm \perp}|^2','interpreter','tex')
% print -djpeg -r300 CavityAbsorptionErez