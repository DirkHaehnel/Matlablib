function fluo2000(action)

%To start the program, simply type: fluo2000
%Author: Joerg Enderlein, 7/1/00 (joerg.enderlein@chemie.uni-regensburg.de) 

global fig h_para h_datawin h_errorwin h_paratext h_plotmodetxt h_plotmode
global brdh brdv hgt top paranum plotmode
global h_offsetchk h_shiftchk h_offsettxt h_shifttxt h_offset h_shift offsetchk shiftchk
global h_dttxt h_dt h_periodtxt h_period h_errtxt h_err dt period err
global offset shift plotdata errdata

brdh = 55;
brdv = 30;
hgt = 25;
top = 35;

if nargin < 1, action = 'initialize'; end;

if strcmp(action,'initialize')
   plotmode = 1;
   paranum = 1;
   offsetchk = 1;
   shiftchk = 1;
   decay = [];
   irf = [];
   dt = 1;
   offset = 0;
   shift = 0;
   period = 200;
   err = 0;
   plotdata = [];
   
   time = 1:128;
   maxtime = 128;
   period = 200;
   deltatau = 1;
   
   set(0,'defaultaxesfontname','times');
   set(0,'defaultaxesfontsize',10);
   set(0,'defaulttextfontname','times');
   set(0,'defaulttextfontsize',10);
   
   fig = figure('Name','FluoFit 2000', 'NumberTitle','off', ...
      'MenuBar', 'none', 'Visible','off', 'BackingStore','off', ...
      'ResizeFcn', 'fluo2000(''resize'')', 'PaperPositionMode', 'auto');
   figure(fig);
   
   h_file = uimenu('Label', '&File');
   uimenu(h_file, 'Label', 'Load IRF', 'Accelerator', 'I', ...
      'Callback', 'fluo2000(''openirf'')');
   uimenu(h_file, 'Label', 'Load decay', 'Accelerator', 'D', ...
      'Callback', 'fluo2000(''opendecay'')');
   uimenu(h_file, 'Label', 'Save', 'Accelerator', 'S', 'Callback', 'fluo2000(''save'')');
   uimenu(h_file, 'Label', 'Print', 'Accelerator', 'P', 'Separator', 'on');
   uimenu(h_file, 'Label', 'Quit', 'Accelerator', 'Q', 'Callback', 'close',...
      'Separator', 'on');
   
   h_fit = uimenu('Label', 'F&it');
   uimenu(h_fit, 'Label', 'Start', 'Accelerator', 'Q', ...
      'Callback', 'fluo2000(''fit'')');
   
   pos = get(fig, 'Position');
   col = get(fig, 'Color');
   h_datawin = subplot('position', [brdh/pos(3) 0.2+2*brdv/pos(4) ...
         (pos(3)-150-brdh)/pos(3) 0.8-2.5*brdv/pos(4)]);
   set(h_datawin, 'Box', 'on', 'xlabel', xlabel('Time (ns)'), 'ylabel', ylabel('Counts'));
   ylabel('Counts');
   h_errorwin = subplot('position', [brdh/pos(3) 0.3*brdv/pos(4) ...
         (pos(3)-150-brdh)/pos(3) 0.2]);
   set(h_errorwin, 'Box', 'on', 'XTick', [], 'ylabel', ylabel('Error'));
   ylabel('Error'); 
   h_plotmodetxt = [uicontrol('Parent', fig, 'Style', 'Text', ...
      'String', 'Lin', 'Position', [pos(3)-125 pos(4)-top-3 25 20], ...
      'HorizontalAlignment','right','BackgroundColor', col),...
      uicontrol('Parent', fig, 'Style', 'Text', ...
      'String', 'Log', 'Position', [pos(3)-55 pos(4)-top-3 25 20], ...
      'HorizontalAlignment','right','BackgroundColor', col)];
   h_plotmode = [uicontrol('Parent', fig, 'Style', 'Radiobutton', ...
      'Position', [pos(3)-95 pos(4)-top 15 20], 'Value', plotmode, ...
      'BackgroundColor', col, 'Callback','fluo2000(''plotlin'')'),...
      uicontrol('Parent', fig, 'Style', 'Radiobutton', ...
      'Position', [pos(3)-25 pos(4)-top 15 20], 'Value', 1-plotmode, ...
      'BackgroundColor', col, 'Callback','fluo2000(''plotlog'')')];
   h_paratext = uicontrol('Parent', fig, 'Style', 'Text', ...
      'String', '# of Parameters', 'Position', [pos(3)-145 pos(4)-top-3-hgt 100 20], ...
      'HorizontalAlignment','right','BackgroundColor', col);
   h_para = uicontrol('Parent', fig, 'Style', 'Edit', ...
      'Position', [pos(3)-35 pos(4)-top-hgt 25 20], 'BackgroundColor', 'w', 'String', 1,...
      'Callback','fluo2000(''paranum'')');
   fluo2000('uicontrol');
end
   
if strcmp(action,'uicontrol')
   pos = get(fig, 'Position');
   col = get(fig, 'Color');
   for j=1:paranum
      h_paratext = [h_paratext uicontrol('Parent', fig, 'Style', 'Text', ...
         'String', ['Tau' num2str(j)], 'Position', [pos(3)-145 pos(4)-top-3-2*j*hgt 75 20], ...
         'HorizontalAlignment','right','BackgroundColor', col)];
      h_paratext = [h_paratext uicontrol('Parent', fig, 'Style', 'Text', ...
         'String', ['Amp' num2str(j)], 'Position', [pos(3)-145 pos(4)-top-3-(2*j+1)*hgt 75 20], ...
         'HorizontalAlignment','right','BackgroundColor', col)];
   end
   h_offsettxt = uicontrol('Parent', fig, 'Style', 'Text', ...
         'String', 'Offset', 'Position', [pos(3)-115 pos(4)-top-3-(2*paranum+2)*hgt 45 20], ...
         'HorizontalAlignment','right','BackgroundColor', col);
   h_shifttxt = uicontrol('Parent', fig, 'Style', 'Text', ...
         'String', 'Col. shift' , 'Position', [pos(3)-115 pos(4)-top-3-(2*paranum+3)*hgt 45 20], ...
         'HorizontalAlignment','right','BackgroundColor', col);
   h_dttxt = uicontrol('Parent', fig, 'Style', 'Text', ...
         'String', 'Channel (ns)', 'Position', [pos(3)-145 pos(4)-top-3-(2*paranum+4.5)*hgt 75 20], ...
         'HorizontalAlignment','right','BackgroundColor', col);
   h_periodtxt = uicontrol('Parent', fig, 'Style', 'Text', ...
         'String', 'Period (ns)' , 'Position', [pos(3)-145 pos(4)-top-3-(2*paranum+5.5)*hgt 75 20], ...
         'HorizontalAlignment','right','BackgroundColor', col);
   h_errtxt = uicontrol('Parent', fig, 'Style', 'Text', ...
         'String', 'Error \chi^2' , 'Position', [pos(3)-145 pos(4)-top-3-(2*paranum+7)*hgt 75 20], ...
         'HorizontalAlignment','right','BackgroundColor', col);
   for j=1:paranum
      h_para = [h_para uicontrol('Parent', fig, 'Style', 'Edit', ...
            'Position', [pos(3)-60 pos(4)-top-2*j*hgt 50 20], ...
            'BackgroundColor', 'w', 'Callback','fluo2000(''paranum'')')];
      h_para = [h_para uicontrol('Parent', fig, 'Style', 'Edit', ...
            'Position', [pos(3)-60 pos(4)-top-(2*j+1)*hgt 50 20], ...
            'BackgroundColor', 'w', 'Callback','fluo2000(''paranum'')')];
   end
   h_offset = uicontrol('Parent', fig, 'Style', 'Edit', ...
      'Position', [pos(3)-60 pos(4)-top-(2*paranum+2)*hgt 50 20], ...
      'BackgroundColor', 'w', 'String', offset,...
      'Callback', 'offset = str2num(get(h_offset,''String''))');
   h_shift = uicontrol('Parent', fig, 'Style', 'Edit', ...
      'Position', [pos(3)-60 pos(4)-top-(2*paranum+3)*hgt 50 20], ...
      'BackgroundColor', 'w', 'String', shift,...
      'Callback', 'shift = str2num(get(h_shift,''String''))');
   h_dt = uicontrol('Parent', fig, 'Style', 'Edit', ...
      'Position', [pos(3)-60 pos(4)-top-(2*paranum+4.5)*hgt 50 20], ...
      'BackgroundColor', 'w', 'String', dt, ...
      'Callback', 'dt = str2num(get(h_dt,''String''))');
   h_period = uicontrol('Parent', fig, 'Style', 'Edit', ...
      'Position', [pos(3)-60 pos(4)-top-(2*paranum+5.5)*hgt 50 20], ...
      'BackgroundColor', 'w', 'String', period, ...
      'Callback', 'period = str2num(get(h_period,''String''))');
   h_err = uicontrol('Parent', fig, 'Style', 'Edit', ...
      'Position', [pos(3)-60 pos(4)-top-(2*paranum+7)*hgt 50 20], ...
      'String', err, 'Enable', 'inactive');
   h_offsetchk = uicontrol('Parent', fig, 'Style', 'Radiobutton', ...
      'Value', offsetchk, 'Position', [pos(3)-135 pos(4)-top-3-(2*paranum+2)*hgt 15 20],...
      'BackgroundColor', col, 'Callback', 'fluo2000(''offsetchk'')');
   h_shiftchk = uicontrol('Parent', fig, 'Style', 'Radiobutton', ...
      'Value', shiftchk, 'Position', [pos(3)-135 pos(4)-top-3-(2*paranum+3)*hgt 15 20],...
      'BackgroundColor', col, 'Callback', 'fluo2000(''shiftchk'')');
end; 	% end initialization routine 

if strcmp(action,'resize')
   pos = get(fig, 'Position');
   set(h_datawin, 'position', [brdh/pos(3) 0.2+2*brdv/pos(4) ...
         (pos(3)-150-brdh)/pos(3) 0.8-2.5*brdv/pos(4)]);
   set(h_errorwin, 'position', [brdh/pos(3) 0.3*brdv/pos(4) ...
         (pos(3)-150-brdh)/pos(3) 0.2]);
   set(h_plotmodetxt(1), 'Position', [pos(3)-125 pos(4)-top-3 25 20]);
   set(h_plotmodetxt(2), 'Position', [pos(3)-55 pos(4)-top-3 25 20]);
   set(h_plotmode(1), 'Position', [pos(3)-95 pos(4)-top 15 20]);
   set(h_plotmode(2), 'Position', [pos(3)-25 pos(4)-top 15 20]);
   set(h_paratext(1), 'Position', [pos(3)-145 pos(4)-top-hgt-3 100 20]);
   set(h_para(1), 'Position', [pos(3)-35 pos(4)-top-hgt 25 20]);
   for j=1:2*paranum
      set(h_paratext(j+1), 'Position', [pos(3)-145 pos(4)-top-3-(j+1)*hgt 75 20]);
      set(h_para(j+1), 'Position', [pos(3)-60 pos(4)-top-(j+1)*hgt 50 20]);
   end
   set(h_offsetchk, 'Position', [pos(3)-135 pos(4)-top-3-(2*paranum+2)*hgt 15 20]);
   set(h_shiftchk, 'Position', [pos(3)-135 pos(4)-top-3-(2*paranum+3)*hgt 15 20]);
   set(h_offsettxt, 'Position', [pos(3)-115 pos(4)-top-3-(2*paranum+2)*hgt 45 20]);
   set(h_shifttxt, 'Position', [pos(3)-115 pos(4)-top-3-(2*paranum+3)*hgt 45 20]);
   set(h_dttxt, 'Position', [pos(3)-145 pos(4)-top-3-(2*paranum+4.5)*hgt 75 20]);
   set(h_periodtxt, 'Position', [pos(3)-145 pos(4)-top-3-(2*paranum+5.5)*hgt 75 20]);
   set(h_errtxt, 'Position', [pos(3)-145 pos(4)-top-3-(2*paranum+7)*hgt 75 20]);
   set(h_offset, 'Position', [pos(3)-60 pos(4)-top-(2*paranum+2)*hgt 50 20]);
   set(h_shift, 'Position', [pos(3)-60 pos(4)-top-(2*paranum+3)*hgt 50 20]);
   set(h_dt, 'Position', [pos(3)-60 pos(4)-top-(2*paranum+4.5)*hgt 50 20]);
   set(h_period, 'Position', [pos(3)-60 pos(4)-top-(2*paranum+5.5)*hgt 50 20]);
   set(h_err, 'Position', [pos(3)-60 pos(4)-top-(2*paranum+7)*hgt 50 20]);
end

if strcmp(action,'openirf')
   [fname, pname] = uigetfile('*.*');
   fin = fopen([pname fname],'r');
   irf = fscanf(fin, '%f', inf);
   fclose(fin);
   plotdata = [irf(:) plotdata];
   time = dt*(1:length(irf));
   subplot(h_datawin);
   tmp = get(h_datawin, 'children');
   set(h_datawin, 'children', plot(time, plotdata, 'erasemode', 'xor'),...
      'xlim', [0 max(time)], 'xlabel', xlabel('Time(ns)'), 'ylabel', ylabel('Counts'));
end

if strcmp(action,'opendecay')
   [fname, pname] = uigetfile('*.*');
   fin = fopen([pname fname],'r');
   tmp = fscanf(fin, '%f', inf);
   fclose(fin);
   plotdata = [plotdata tmp(:)];
   time = dt*(1:length(plotdata(:,1)));
   subplot(h_datawin);
   set(h_datawin, 'children', plot(time, plotdata, 'erasemode', 'xor' ),...
      'xlim', [0 max(time)], 'xlabel', xlabel('Time(ns)'), 'ylabel', ylabel('Counts'));
end

if strcmp(action,'save')
   [fname, pname] = uiputfile('*.*');
   disp(fname);
end

if strcmp(action,'plotlin')
   plotmode = 1;
   set(h_plotmode(1),'Value',plotmode);
   set(h_plotmode(2),'Value',1-plotmode);
   set(h_datawin, 'yscale', 'linear');
end

if strcmp(action,'plotlog')
   plotmode = 0;
   set(h_plotmode(1),'Value',plotmode);
   set(h_plotmode(2),'Value',1-plotmode);
   set(h_datawin, 'yscale', 'log');
end

if strcmp(action,'paranum')
   num = round(str2num(get(h_para(1),'String')));
   if num>0 & num<5 & ~(num==paranum)
      for j=1:paranum
         delete(h_paratext(2*j));
         delete(h_paratext(2*j+1));
         delete(h_para(2*j));
         delete(h_para(2*j+1));
      end
      h_paratext(2:end) = [];
      h_para(2:end) = [];
      delete(h_offsetchk);
      delete(h_shiftchk);
      delete(h_offsettxt);
      delete(h_shifttxt);
      delete(h_dttxt);
      delete(h_periodtxt);
      delete(h_errtxt);
      delete(h_offset);
      delete(h_shift);
      delete(h_dt);
      delete(h_period);
      delete(h_err);
      paranum = num;
      set(h_para(1), 'String', paranum);
      fluo2000('uicontrol');
   end
end

if strcmp(action,'offsetchk')
   offsetchk = get(h_offsetchk, 'Value');
   if offsetchk==0
      col = get(fig, 'Color');
      set(h_offset, 'BackgroundColor', col, 'String', 0);
   else
      set(h_offset, 'BackgroundColor', 'w');
   end
end

if strcmp(action,'shiftchk')
   shiftchk = get(h_shiftchk, 'Value');
   if shiftchk==0
      col = get(fig, 'Color');
      set(h_shift, 'BackgroundColor', col, 'String', 0);
   else
      set(h_shift, 'BackgroundColor', 'w');
   end
end

if strcmp(action,'fit')
   subplot(h_datawin);
   y = decay(:,1)';
   irf = irf(:)';
   init = 1;
   m = paranum;
   p = period;
   
   p = p/dt;
   n = length(irf); 
   tp = 1:p;
   t = 1:length(y);
   
   if shiftchk
      sh_min = -5;
      sh_max = 5;
   else
      sh_min = 0;
      sh_max = 0;
   end
   if init
      param(1) = m;
      for i=2:m 
         param(i) = -param(i-1)*(m-i+1)/i;
      end;
      B = y(m+1:n)';
      for i=1:m
         B = B + param(i)^2*y(m+1-i:n-i)';
      end
      for i=1:m
         tmp = -param(i)*y(m+1:n-i)';
         for j=1:m-i
            tmp = tmp + param(j)*param(j+i)*y(m+1-j:n-j-i)';
         end
         B = [[tmp; zeros(i,1)] B [zeros(i,1); tmp]];
      end
      d = -m:m;
      V = spdiags(B, d, n-m, n-m);
      D = [y(m:n-1)'];
      for i=2:m
         D = [D y(m+1-i:n-i)'];
      end
      D = D + 1e-10*(D==0);
      D(:,m+1:2*m+1)=ones(n-m,m+1);
      yy = y(m+1:n)';
      for k=sh_min:sh_max
         irs=irf(rem(rem(t-k-1, n)+n,n)+1);
         for i=1:m
            D(:,m+i) = irs(m+2-i:n+1-i)';
         end
         VD = V\D;
         param = (VD'*D)\(VD'*yy);
         tau = roots([param(m:-1:1)' -1]);
         if ~isempty(tau)
            for i=1:m
               tmp = (1-tau(i)./tau);
               tmp(i) = [];
               A(i) = polyval(param(2*m:-1:m+1), tau(i))/prod(tmp);
            end
            tau = 1./log(tau);
            offset = param(2*m+1)/(1-sum(param(1:m)));
            z = A*diag(1./(1-exp(-p./tau)))*exp(-(1./tau)*(tp-1));
            z = convol(irs, z) + offset;
            tmp = ((y-z)./abs(z))*(y-z)'/(n-m);
            chi(k-sh_min+1) = tmp;
         else
            chi(k-sh_min+1) = nan;
         end
      end
      [chi, c] = min(chi);
      c = c + sh_min - 1;
      irs=irf(rem(rem(t-c-1, n)+n,n)+1);
      for i=1:m
         D(:,m+i) = irs(m+2-i:n+1-i)';
      end
      VD = V\D;
      param = (VD'*D)\(VD'*yy);
      tau = roots([param(m:-1:1)' -1]);
      for i=1:m
         tmp = (1-tau(i)./tau);
         tmp(i) = [];
         A(i) = polyval(param(2*m:-1:m+1), tau(i))/prod(tmp);
      end
      tau = abs(1./log(tau));
      offset = param(2*m+1)/(1-sum(param(1:m)));
      z = A*diag(1./(1-exp(-p./tau)))*exp(-(1./tau)*(tp-1));
      z = convol(irs, z) + offset;
   else
      tau = tau(:)
      m = length(tau);
      x = diag(1./(1-exp(-p./tau)))*exp(-(1./tau)*(tp-1));
      z = zeros(m, n);
      for j=1:m
         z(j,:) = convol(irf, x(j,:));
      end
      A = y/z;
      z = A*z;
      c = 0;
      offset = options(2);
   end
   if offsetchk offset=abs(offset); else offset = 0; end
   tmp = get(h_datawin,'children');
   subplot(h_datawin);
   plot(t,irf,'g',t,y,'r',t,z,'erasemode','xor');
   xlabel('Time (ns)'); ylabel('Counts');
   param = [c; offset; tau];
   % Decay times and Offset are assumed to be positive.
   paramin = [-Inf, 0, zeros(1,length(tau))];
   paramax = [Inf, Inf, Inf*ones(1,length(tau))];
   if ~shiftchk  paramin(1) = 0; paramax(1) = 0; end
   if ~offsetchk  paramax(2) = 0; end
   [param, dparam] = simplex('fluofitfun', param, paramin, paramax, 1e-8, [], irf(:), y(:), p);
   c = param(1);
   dc = dparam(1);
   offset = param(2);
   doffset = dparam(2);
   tau = param(3:length(param));
   dtau = dparam(3:length(param));
   x = diag(1./(1-exp(-p./tau)))*exp(-(1./tau)*(tp-1));
   irs = (1-c+floor(c))*irf(rem(rem(t-floor(c)-1, n)+n,n)+1) + (c-floor(c))*irf(rem(rem(t-ceil(c)-1, n)+n,n)+1);
   for j =1:m
      z(j,:) = convol(irs, x(j,:));
   end
   A = (z*z')\(z*(y'-offset));
   z = A'*z + offset;
   chi = ((y-z)./abs(z))*(y-z)'/(n-m);
   t = dt*t;
   tau = dt*tau;
   c = dt*c;
   dtau = dt*dtau;
   dc = dt*dc;
   set(plotdata, 'ydata', z);
   v = axis;
   v(1) = min(t);
   v(2) = max(t);
   axis(v);
   xlabel('Time in ns');
   ylabel('Log Count');
   subplot(h_errorwin);
   plot(t,(y-z)./sqrt(abs(z)));
   v = axis;
   v(1) = min(t);
   v(2) = max(t);
   axis(v);
   xlabel('Time in ns');
   ylabel('Residue');
end
