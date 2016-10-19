function scan2000(action, varin)

%To start the program, simply type: scan2000
%Author: Joerg Enderlein, 11/11/00 (joerg.enderlein@chemie.uni-regensburg.de) 

%handles:
global fig h_intwin h_tauwin h_intbar h_taubar h_intim h_tauim 
global h_para h_paratext h_plotmodetxt h_plotmode h_button
%window border values:
global brdh brdv height txtwidth editwidth split axisbrd
%data:
global path para buttonname plotmode x y tcspc tag

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
   
   paraname = {'Bin width', 'Offset', 'X-pixels', 'Y-pixels', 'min. Int.', 'max. Int.', 'min. Tau', 'max. Tau', ...
         'min. Cut', 'max. Cut'};
   para = [5000 7 500 500 0 0 0 0 0 0];
	buttonname = {'Reevaluate', 'Lie filter', 'Molecules'};
	buttonfunc = {'reeval', 'lie', 'locate'};
     
   plotmode = 1;   
   
	scrsz = get(0,'ScreenSize');   
   fig = figure('Name','Scan 2000', 'NumberTitle','off', ...
      'MenuBar', 'none', 'Visible','off', 'BackingStore','off', 'Position', [scrsz(3)/10 scrsz(4)/10 4*scrsz(3)/5 4*scrsz(4)/5], ...
      'ResizeFcn', 'scan2000(''resize'')', 'PaperPositionMode', 'auto');
   figure(fig);
   
   %File menu
   h_file = uimenu('Label', '&File');
   uimenu(h_file, 'Label', 'Load Data', 'Accelerator', 'D', ...
      'Callback', 'scan2000(''open'')');
   uimenu(h_file, 'Label', 'Save', 'Accelerator', 'S', 'Callback', 'scan2000(''save'')');
   uimenu(h_file, 'Label', 'Print', 'Accelerator', 'P', 'Separator', 'on');
   uimenu(h_file, 'Label', 'Quit', 'Accelerator', 'Q', 'Callback', 'close',...
      'Separator', 'on');
   
   %Image menu
   h_image = uimenu('Label', '&Image');
   uimenu(h_image, 'Label', 'Zoom', 'Accelerator', 'Z', ...
      'Callback', 'scan2000(''zoom'')');
  
   %Data menu
   h_data = uimenu('Label', '&Data');
   uimenu(h_data, 'Label', 'Lie filter', 'Accelerator', 'L', ...
      'Callback', 'scan2000(''lie'')');
   uimenu(h_data, 'Label', 'Locate SM', 'Accelerator', 'S', ...
      'Callback', 'scan2000(''locate'')');
     
   %Help menu
   h_help = uimenu('Label', '&Help');
   uimenu(h_help, 'Label', 'Help', 'Accelerator', 'L', ...
      'Callback', 'scan2000(''help'')');
   uimenu(h_help, 'Label', 'About', 'Accelerator', 'A', ...
      'Callback', 'scan2000(''about'')', 'Separator', 'on');
   
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
         'BackgroundColor', 'w', 'String', para(j), 'Callback','scan2000(''para'')')];
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
         'Value', plotmode, 'BackgroundColor', col, 'Callback','scan2000(''plotlin'')'),...
         uicontrol('Parent', fig, 'Style', 'Radiobutton', 'Position', ...
         [pos(3)-brdh-editwidth pos(4)-brdh-3-(j+1)*height-j*split editwidth height], ...
         'Value', 1-plotmode, 'BackgroundColor', col, 'Callback','scan2000(''plotlog'')')];
   k = length(para)+2;
   h_button = [];
   for j = 1:length(buttonname)
      h_button = [h_button, uicontrol('Parent', fig, 'Style', 'pushbutton', 'Position', ...
            [pos(3)-brdh-editwidth-textwidth pos(4)-brdh-3-(k+j)*height-(k+j-1)*split editwidth+textwidth height], ...
            'String', buttonname{j}, 'BackgroundColor', col, 'Callback',['scan2000('''  buttonfunc{j} ''')'])]; 
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
      fin = fopen([pname fname],'r');
      if (fin==-1)
         errordlg('Cannot open speciefied file. Please try again.');
      else        
         h = waitbar(0,'Please wait ...'); drawnow
         fseek(fin, 48+256, 0);
         NChannels = fread(fin, 1, 'int32');
         NCurves = fread(fin, 1, 'int32');
         ActiveCurve = fread(fin, 1, 'int32');
         MeasMode = fread(fin, 1, 'int32');
         HistMode = fread(fin, 1, 'int32');
         Range = fread(fin, 1, 'int32');
         Offset = fread(fin, 1, 'int32');
         AcquTime = fread(fin, 1, 'int32');
         LinLog = fread(fin, 1, 'int32');
         MinAx = fread(fin, 1, 'int32');
         MaxAx = fread(fin, 1, 'int32');
         MinAxCnt = fread(fin, 1, 'int32');
         MaxAxCnt = fread(fin, 1, 'int32');
         CFD0 = fread(fin, 1, 'int32');
         CFDmin = fread(fin, 1, 'int32');
         Sync = fread(fin, 1, 'int32');
         Resolution = fread(fin, 1, 'float');
         GlobClock = fread(fin, 1, 'int32');
         NCounts = fread(fin, 1, 'uint32');
         y = fread(fin, NCounts, 'uint32');   
         fclose(fin);
         waitbar(0.25);
         flag = floor(y/2^30);
         y = y - floor(y/2^28)*2^28;
         x = floor(y/2^16);
         y = y - x*2^16;
         x = 2^12 + 1 - x; 
         y = y + 2^16*cumsum(~(flag==1));
         y(~(flag==1))=[];
         x(~(flag==1))=[];
         y = y - y(1);
         y = tttr2bin(y,para(1));
         waitbar(0.5);
         t = reshape(1:para(3)*para(4),para(3),para(4));
         t(:,2:2:end) = t(end:-1:1,2:2:end);
         t = reshape(t,1,para(3)*para(4));
         j = 0;
         tmp = y(j+t(1:para(3)*(para(4)-1)))'*y(j+t((para(3)+1):para(3)*para(4)));
         for j=1:20
            tmp = [tmp, y(j+t(1:para(3)*(para(4)-1)))'*y(j+t((para(3)+1):para(3)*para(4)))]; 
            if tmp(end)<tmp(end-1) break; end
         end
         [tmp, tmp] = max(tmp);      
         para(2) = tmp-1;
         waitbar(0.75);
         tmp = cumsum(y);
         x = cumsum([0; x]);
         x = (x(tmp+1)-x([1; tmp(1:end-1)]))./(y+(y==0));
         x = x-min(x);
         tag = reshape(y(para(2)+(1:para(3)*para(4))),para(3),para(4));
         tcspc = reshape(x(para(2)+(1:para(3)*para(4))),para(3),para(4));   
         tag(:,2:2:end) = tag(end:-1:1,2:2:end);
         tcspc(:,2:2:end) = tcspc(end:-1:1,2:2:end);   
         para(5) = min(min(tag));
         para(6) = max(max(tag));
         para(7) = min(min(tcspc));
         para(8) = max(max(tcspc));
         para(9) = para(5);
         para(10) = para(6);
         for j=1:length(para)
            set(h_para(j),'String',num2str(para(j)));
         end
         close(h);
         scan2000('reeval');
      end
   end
end

if strcmp(action,'reeval')
   if ~isempty(x)
      subplot(h_intwin);
      plt = tag; plt(plt<para(5)) = para(5); plt(plt>para(6)) = para(6);
      imagesc(plt); axis image
      set(h_intwin, 'ButtonDownFcn', 'scan2000(''zoom'')');
      subplot(h_intbar);
      tmp = para(6)-para(5);
      tmp = para(5) + (0:100)/100*tmp;
      imagesc(tmp,[0 0.05*(para(6)-para(5))],[tmp; tmp]); axis image;
      set(h_intbar,'ytick',[]);
      subplot(h_tauwin);
      plt = tcspc; 
      plt(plt<para(7)) = para(7); plt(plt>para(8)) = para(8);
      plt(tag<para(9) & tag<=para(10)) = para(7);
      imagesc(plt); axis image
      set(h_tauwin,'ytick',[]);
      set(h_tauwin, 'ButtonDownFcn', 'scan2000(''zoom'')');
      subplot(h_taubar);
      tmp = para(8)-para(7);
      tmp = para(7) + (0:100)/100*tmp;
      imagesc(tmp,[0 0.05*(para(8)-para(7))],[tmp; tmp]); axis image;
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
   tmp = round(min(min(tcspc(tag>=para(9) & tag<=para(10)))));
   if tmp>para(7) 
      para(7) = tmp; 
      set(h_para(7),'String',num2str(para(7)));
   end
   tmp = round(max(max(tcspc(tag>=para(9) & tag<=para(10)))));
   if tmp<para(8) 
      para(8) = tmp; 
      set(h_para(8),'String',num2str(para(8)));
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
