function scan2001(action, varin)

%To start the program, simply type: scan2001
%Author: Joerg Enderlein, 11/11/00 (joerg.enderlein@chemie.uni-regensburg.de) 

%handles:
global fig h_intwin h_tauwin h_intbar h_taubar h_intim h_tauim 
global h_para h_paratext h_plotmodetxt h_plotmode h_button
%window border values:
global brdh brdv height txtwidth editwidth split axisbrd
%data:
global path para buttonname plotmode x tcspc0 tag0 tcspc tag

brdh = 5; % horizontal border
brdv = 5; % vertical border
axisbrd = 30; % border for axis labels
textwidth = 50; % width of text
editwidth = 40; % width of edit window
height = 25; %height of editwindow
split = 5; % vertical distance between editwindows

if nargin < 1, action = 'initialize'; end;

if strcmp(action,'initialize')
   set(0,'defaultaxesfontname','times');
   set(0,'defaultaxesfontsize',12);
   set(0,'defaulttextfontname','times');
   set(0,'defaulttextfontsize',12);
   
   tag0 = [];
   tcspc0 = [];   
   tag = [];   
   tcspc = [];   
   x = [];
   
   paraname = {'Bin width', 'X-pixels', 'Y-pixels', 'min. Int.', 'max. Int.', 'min. Tau', 'max. Tau', ...
         'min. Cut', 'max. Cut'};
   para = [5000 500 500 0 0 0 0 0 0];
	buttonname = {'Reevaluate', 'Lie filter', 'Molecules'};
	buttonfunc = {'reeval', 'lie', 'locate'};
     
   plotmode = 1;   
   
	scrsz = get(0,'ScreenSize');   
   fig = figure('Name','Scan 2001', 'NumberTitle','off', ...
      'MenuBar', 'none', 'Visible','off', 'BackingStore','off', 'Position', [scrsz(3)/10 scrsz(4)/10 4*scrsz(3)/5 4*scrsz(4)/5], ...
      'ResizeFcn', 'scan2001(''resize'')', 'PaperPositionMode', 'auto');
   figure(fig);
   
   %File menu
   h_file = uimenu('Label', '&File');
   uimenu(h_file, 'Label', 'Load Data', 'Accelerator', 'D', ...
      'Callback', 'scan2001(''open'')');
   uimenu(h_file, 'Label', 'Save', 'Accelerator', 'S', 'Callback', 'scan2001(''save'')');
   uimenu(h_file, 'Label', 'Print', 'Accelerator', 'P', 'Separator', 'on');
   uimenu(h_file, 'Label', 'Quit', 'Accelerator', 'Q', 'Callback', 'close',...
      'Separator', 'on');
   
   %Image menu
   h_image = uimenu('Label', '&Image');
   uimenu(h_image, 'Label', 'Zoom', 'Accelerator', 'Z', ...
      'Callback', 'scan2001(''zoom'')');
  
   %Data menu
   h_data = uimenu('Label', '&Data');
   uimenu(h_data, 'Label', 'Lie filter', 'Accelerator', 'L', ...
      'Callback', 'scan2001(''lie'')');
   uimenu(h_data, 'Label', 'Locate SM', 'Accelerator', 'S', ...
      'Callback', 'scan2001(''locate'')');
     
   %Help menu
   h_help = uimenu('Label', '&Help');
   uimenu(h_help, 'Label', 'Help', 'Accelerator', 'L', ...
      'Callback', 'scan2001(''help'')');
   uimenu(h_help, 'Label', 'About', 'Accelerator', 'A', ...
      'Callback', 'scan2001(''about'')', 'Separator', 'on');
   
   pos = get(fig, 'Position');
   col = get(fig, 'Color');
   
   rngv = pos(4)-3*brdv-2*axisbrd;
   rngh = pos(3)-5*brdh-textwidth-editwidth-split-axisbrd;
   h_intwin = subplot('position', [(axisbrd+brdh)/pos(3) (2*brdv+2*axisbrd+0.1*rngv)/pos(4) rngh/2/pos(3) 0.9*rngv/pos(4)]);
   set(h_intwin, 'Box', 'on', 'xtick', [], 'ytick', []);
   h_tauwin = subplot('position', [(axisbrd+2*brdh+rngh/2)/pos(3) (2*brdv+2*axisbrd+0.1*rngv)/pos(4) rngh/2/pos(3) 0.9*rngv/pos(4)]);
   set(h_tauwin, 'Box', 'on', 'xtick', [], 'ytick', []);
   h_intbar = subplot('position', [(axisbrd+brdh)/pos(3) (axisbrd+brdv)/pos(4) rngh/2/pos(3) 0.1*rngv/pos(4)]);
   set(h_intbar, 'Box', 'on', 'xtick', [], 'ytick', []);
   h_taubar = subplot('position', [(axisbrd+2*brdh+rngh/2)/pos(3) (axisbrd+brdv)/pos(4) rngh/2/pos(3) 0.1*rngv/pos(4)]);
   set(h_taubar, 'Box', 'on', 'xtick', [], 'ytick', []);
   h_paratext = [];
   h_para = [];
   for j=1:length(para)
	   h_paratext = [h_paratext, uicontrol('Parent', fig, 'Style', 'Text', 'String', paraname{j}, 'Position', ...
         [pos(3)-brdh-textwidth-editwidth-split pos(4)-brdh-8-j*height-(j-1)*split textwidth height], ...
         'HorizontalAlignment','right','BackgroundColor', col)];
	   h_para = [h_para, uicontrol('Parent', fig, 'Style', 'Edit', 'Position', ...
   	   [pos(3)-brdh-editwidth pos(4)-brdh-3-j*height-(j-1)*split editwidth height], ...
         'BackgroundColor', 'w', 'String', para(j), 'Callback','scan2001(''para'')')];
	end
   j = length(para)+1;
   h_plotmodetxt = [uicontrol('Parent', fig, 'Style', 'Text', 'String', 'Lin', 'Position', ...
         [pos(3)-brdh-textwidth-editwidth-split pos(4)-brdh-8-j*height-(j-1)*split textwidth height], ...
         'HorizontalAlignment','right','BackgroundColor', col),...
         uicontrol('Parent', fig, 'Style', 'Text', 'String', 'Log', 'Position', ...
         [pos(3)-brdh-textwidth-editwidth-split pos(4)-brdh-8-(j+1)*height-j*split textwidth height], ...
         'HorizontalAlignment','right','BackgroundColor', col)];
   h_plotmode = [uicontrol('Parent', fig, 'Style', 'Radiobutton', 'Position', ...
         [pos(3)-brdh-editwidth pos(4)-brdh-3-j*height-(j-1)*split editwidth height], ...
         'Value', plotmode, 'BackgroundColor', col, 'Callback','scan2001(''plotlin'')'),...
         uicontrol('Parent', fig, 'Style', 'Radiobutton', 'Position', ...
         [pos(3)-brdh-editwidth pos(4)-brdh-3-(j+1)*height-j*split editwidth height], ...
         'Value', 1-plotmode, 'BackgroundColor', col, 'Callback','scan2001(''plotlog'')')];
   k = length(para)+2;
   h_button = [];
   for j = 1:length(buttonname)
      h_button = [h_button, uicontrol('Parent', fig, 'Style', 'pushbutton', 'Position', ...
            [pos(3)-brdh-editwidth-textwidth pos(4)-brdh-3-(k+j)*height-(k+j-1)*split editwidth+textwidth height], ...
            'String', buttonname{j}, 'BackgroundColor', col, 'Callback',['scan2001('''  buttonfunc{j} ''')'])]; 
   end
end

if strcmp(action,'resize')
   pos = get(fig, 'position');
   rngv = pos(4)-3*brdv-2*axisbrd;
   rngh = pos(3)-5*brdh-textwidth-editwidth-split-axisbrd;
   set(h_intwin, 'position', [(axisbrd+brdh)/pos(3) (2*brdv+2*axisbrd+0.1*rngv)/pos(4) rngh/2/pos(3) 0.9*rngv/pos(4)]);
   set(h_tauwin, 'position', [(axisbrd+2*brdh+rngh/2)/pos(3) (2*brdv+2*axisbrd+0.1*rngv)/pos(4) rngh/2/pos(3) 0.9*rngv/pos(4)]);
   set(h_intbar, 'position', [(axisbrd+brdh)/pos(3) (axisbrd+brdv)/pos(4) rngh/2/pos(3) 0.1*rngv/pos(4)]);
   set(h_taubar, 'position', [(axisbrd+2*brdh+rngh/2)/pos(3) (axisbrd+brdv)/pos(4) rngh/2/pos(3) 0.1*rngv/pos(4)]);
   for j=1:length(para)
      set(h_paratext(j), 'Position', [pos(3)-2*brdh-textwidth-editwidth-split pos(4)-brdh-8-j*height-(j-1)*split textwidth height]);
      set(h_para(j), 'Position', [pos(3)-2*brdh-editwidth pos(4)-brdh-8-j*height-(j-1)*split editwidth height]);
   end
   for k = 1:2
      j = length(para)+k;
      set(h_plotmodetxt(k), 'Position', ...
         [pos(3)-brdh-textwidth-editwidth-split pos(4)-brdh-8-j*height-(j-1)*split textwidth height]);
      set(h_plotmode(k), 'Position', ...
         [pos(3)-brdh-editwidth pos(4)-brdh-3-j*height-(j-1)*split editwidth height]);
   end
   k = length(para)+2;
   for j = 1:length(buttonname)
      set(h_button(j), 'Position', ...
            [pos(3)-brdh-editwidth-textwidth pos(4)-brdh-3-(k+j)*height-(k+j-1)*split editwidth+textwidth height]); 
   end
end

if strcmp(action,'open')
   [fname, pname] = uigetfile([path '*.t3r']);
   if ~(fname==0)
      path = pname;
      [tag0, tcspc0, x] = ScanRead([pname fname], para(1:3));
      if ~isempty(tag)   
         tag = tag0;
         tcspc = tcspc0;
         para(4) = min(min(tag));
         para(5) = max(max(tag));
         para(6) = min(min(tcspc));
         para(7) = max(max(tcspc));
         para(8) = para(4);
         para(9) = para(5);
         for j=1:length(para)
            set(h_para(j),'String',num2str(para(j)));
         end
         scan2001('reeval');
      end
   end
end

if strcmp(action,'reeval')
   if ~isempty(tag)
      subplot(h_intwin);
      plt = tag; plt(plt<para(4)) = para(4); plt(plt>para(5)) = para(5);
      imagesc(plt); axis image
      set(h_intwin, 'ButtonDownFcn', 'scan2001(''zoom'')');
      subplot(h_intbar);
      tmp = para(5)-para(4);
      tmp = para(4) + (0:100)/100*tmp;
      imagesc(tmp,[0 0.05*(para(5)-para(4))],[tmp; tmp]); axis image;
      set(h_intbar,'ytick',[]);
      subplot(h_tauwin);
      plt = tcspc; 
      plt(plt<para(6)) = para(6); plt(plt>para(7)) = para(7);
      plt(tag<para(8) & tag<=para(9)) = para(6);
      imagesc(plt); axis image
      set(h_tauwin,'ytick',[]);
      set(h_tauwin, 'ButtonDownFcn', 'scan2001(''zoom'')');
      subplot(h_taubar);
      tmp = para(7)-para(6);
      tmp = para(6) + (0:100)/100*tmp;
      imagesc(tmp,[0 0.05*(para(7)-para(6))],[tmp; tmp]); axis image;
      set(h_taubar,'ytick',[]);
   end      
end

if strcmp(action,'plotlin')
   if plotmode==0
      plotmode = 1;
      set(h_plotmode(1),'Value',plotmode);
      set(h_plotmode(2),'Value',1-plotmode);
      subplot(h_tauwin);
      imagesc(tcspc); axis image; 
      set(h_tauwin,'ytick',[]);
      subplot(h_taubar);
      tmp = max(max(tcspc(isfinite(tcspc))))-min(min(tcspc(isfinite(tcspc))));
      tmp = min(min(tcspc(isfinite(tcspc)))) + (0:100)/100*tmp;
      imagesc(tmp,[0 0.08*tmp(end)],[tmp; tmp]); axis image;
      set(h_taubar,'ytick',[]);
   end
end

if strcmp(action,'plotlog')
   if plotmode==1
      plotmode = 0;
      set(h_plotmode(1),'Value',plotmode);
      set(h_plotmode(2),'Value',1-plotmode);
      subplot(h_tauwin);
      plt = log(tcspc);
      plt(~isfinite(plt)) = 0;
      imagesc(plt); axis image; 
      set(h_tauwin,'ytick',[]);
      subplot(h_taubar);
      tmp = max(max(plt))-min(min(plt));
      tmp = min(min(plt)) + (0:100)/100*tmp;
      imagesc(tmp,[0 0.08*tmp(end)],[tmp; tmp]); axis image;
      set(h_taubar,'ytick',[]);
   end
end

if strcmp(action,'para')
   for j=1:length(para)
      para(j) = str2num(get(h_para(j),'String'));
   end
   tmp = round(min(min(tcspc(tag>=para(8) & tag<=para(9)))));
   if tmp>para(6) 
      para(6) = tmp; 
      set(h_para(6),'String',num2str(para(6)));
   end
   tmp = round(max(max(tcspc(tag>=para(8) & tag<=para(9)))));
   if tmp<para(7) 
      para(7) = tmp; 
      set(h_para(7),'String',num2str(para(7)));
   end
end

if strcmp(action,'save')
   [fname, pname] = uiputfile('*.*');
   disp(fname);
end

if strcmp(action,'zoom')
   h = gcbo;
   waitforbuttonpress
   cp1 = get(h,'CurrentPoint');
   rbbox
   cp2 = get(h,'CurrentPoint');
   cp1 = cp1(1,1:2);
   cp2 = cp2(1,1:2);   
	xl = get(h,'xlim');
	yl = get(h,'ylim');     
   if xl(1)<cp1(1) & cp1(1)<xl(2) & yl(1)<cp1(2) & cp1(2)<yl(2)
      if sum(abs(cp2-cp1))<1
         subplot(h_intwin);
         axis([0 para(3) 0 para(4)]);   
         subplot(h_tauwin);
         axis([0 para(3) 0 para(4)]);   
      else
         p = [cp1(1) cp2(1) cp1(2) cp2(2)];
         if p(1)>p(2) p(1:2)=p(2:-1:1); end
         if p(3)>p(4) p(3:4)=p(4:-1:3); end
         if p(1)<0 p(1)=0; end
         if p(2)>para(3) p(2)=para(3); end
         if p(3)<0 p(3)=0; end
         if p(4)>para(4) p(4)=para(4); end
         subplot(h_intwin);
         axis(p);
         subplot(h_tauwin);
         axis(p);
      end
   end
end

