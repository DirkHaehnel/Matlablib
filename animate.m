function animate( cmd )

% animate( stack )
%
% animated display of the movie in a window, with minimal user controls.
% Click and drag the mouse over the image to control the display:
%    horizontal motion changes contrast 
%    vertical motion changes brightness
% example: click and drag the mouse left to increase contrast
%
% Copyright EMBL, 2002-2004
%
% Please, help us improve this software with feedback/bugs/suggestions !
% This software is provided at no cost by a public research institution. 
% However, postcards are always welcome!
%
% Francois Nedelec 
% Cell Biology and Biophysics, EMBL; Meyerhofstrasse 1; 69117 Heidelberg; Germany
% http://www.embl.org
% http://www.cytosim.org

%===================== handle callbacks
persistent MD;

if isstr( cmd )
        
    switch cmd
            
        case 'play'
            
            set( MD.hPlay, 'String', 'Stop', 'Callback', [mfilename ' stop'] );
            
            if ( ~ MD.play )
                MD.play = true;
            
                % rewind the movie when at end of it
                if MD.frame >= MD.nbFrames
                    MD.frame = 0;
                end

                %start the timer:
                start(MD.timer);
            end
                        
        case 'next'
            
            if MD.frame >= MD.nbFrames
                if MD.loop
                    MD.frame = 0;
                else
                    animate('stop');
                    return;
                end
            end
            
            MD.frame = MD.frame + 1;
            showImage(MD);
            return;
    
        case 'stop'
            
            stop(MD.timer);
            set( MD.hPlay, 'String', 'Play', 'Callback', [mfilename ' play'] );
            MD.play = false;
            
            
        case 'loop'
            
            MD.loop = ~ MD.loop;
            
        case 'faster'
            
            if ( MD.period > 0.005 )
                MD.period = MD.period / 2;
                stop( MD.timer );
                set( MD.timer, 'Period', MD.period );
                if ( MD.play )
                    start( MD.timer )
                end
            end
            if ( ~ MD.play )
                animate('play');
            end
            
        case 'frameslider'
            MD.frame = round( get( MD.hFrame, 'Value' ) );
            showImage(MD);
            
        case 'Cmin'
            MD.Cmin  = round(get( MD.hCmin, 'Value'));
            set( MD.hCminText, 'String', num2str(MD.Cmin));
            if ( MD.Cmax > MD.Cmin )
                MD.Cscale = 256 / ( MD.Cmax - MD.Cmin );
            end
            showImage(MD);

        case 'Cmax'
            MD.Cmax  = round(get( MD.hCmax, 'Value'));
            set( MD.hCmaxText, 'String', num2str(MD.Cmax));
            if ( MD.Cmax > MD.Cmin )
                MD.Cscale = 256 / ( MD.Cmax - MD.Cmin );
            end
            showImage(MD);
                 
        case 'mouse_down'
            
            MD.mouse = get(0, 'pointerlocation');
            MD.Cmin2 = MD.Cmin;
            MD.Cmax2 = MD.Cmax;
            MD.mode  = 1;
           
        case 'mouse_up'
        
            MD.mode = 0;
            
        case 'mouse_motion'
                        
            if ( MD.mode == 1 )          % button down: modify contrast / brightness
                               
                pos.fig = get( gcbf, 'position');
                pos.mouse = get(0, 'pointerlocation');
                pos.relative = ( pos.mouse(1:2) - MD.mouse(1:2) )./ pos.fig(3:4);
                
                dc = pos.relative * ( MD.Cmax2 - MD.Cmin2 );
                MD.Cmin = round( MD.Cmin2 - dc(1) + dc(2) );
                MD.Cmax = round( MD.Cmax2 + dc(1) + dc(2) );
                
                set( MD.hCmin, 'Value', MD.Cmin);
                set( MD.hCminText, 'String', num2str(MD.Cmin));

                set( MD.hCmax, 'Value', MD.Cmax);
                set( MD.hCmaxText, 'String', num2str(MD.Cmax));
                if ( MD.Cmax > MD.Cmin )
                    MD.Cscale = 256 / ( MD.Cmax - MD.Cmin );
                end

                showImage(MD);
            end
            
        case 'key_pressed'
        
            if get(gcf, 'CurrentCharacter') == 'q'
                stop(MD.timer);
                close( gcf );
            end
            
        otherwise
            error('command unknown');
    end
    return;
end

    
%===================== Variables & Window setup:

if ~ isfield( cmd, 'data' )     
    error('missing .data field')
end

if length( cmd ) < 2     
    error('not a movie file')
end

if isfield( cmd, 'filename' )
    [pathstr, name, ext,versn] = fileparts(cmd(1).filename);
    figname = name;
else
    figname = inputname( 1 );
end

%==============  GLOBAL setup

MD.play      = false;
MD.loop      = false;
MD.frame     = 1;
MD.movie     = cmd;
MD.nbFrames  = length( cmd );
MD.period    = 0.25;

MD.mode      = 0;
MD.mouse     = [ 0 , 0 ];

MD.timer = timer('TimerFcn', [mfilename ' next'], 'Period', MD.period, 'ExecutionMode', 'fixedDelay' );

%================    perform a basic auto-scale if needed:

if ( nargin < 2 ) | isempty( minmax )
    for i=1:length(MD.movie)
        im = MD.movie(i).data;
        mn(i) = double( min(min( im )) );
        mx(i) = double( max(max( im )) );
    end
    MD.Cmax = max( mx );
    MD.Cmin = min( mn );
else           % scaling is given as an argument:
    MD.Cmin = minmax(1);
    MD.Cmax = minmax(2);
end

if ( MD.Cmax > MD.Cmin )
    MD.Cscale = 256 / ( MD.Cmax - MD.Cmin );
else
    MD.Cscale = 1;
end

%================   calculate figure size:

imsize   = size( cmd(1).data );
scrsz    = get(0,'ScreenSize');
scrsz    = scrsz(3:4) - [ 10 20 ];
FigScale = min( scrsz(2:-1:1)  ./ imsize );
        
if ( FigScale > 1 )  FigScale = floor( FigScale ); end
if ( min(imsize) > 100 ) & ( FigScale > 1 ) FigScale = 1; end
if ( FigScale > 4 )  FigScale = 4; end
if ( FigScale < 1 )  FigScale = 1 / ceil( 1/FigScale ); end
figsize = FigScale * imsize(2:-1:1);
if ( figsize(1) < 420 ) figsize(1) = 420; end

%================    position the figure:

figpos     = [scrsz(1)/2+10, (scrsz(2) - figsize(2))/2, figsize + [ 0 50 ] ];
hFig = figure( 'Position', figpos, 'Name', figname, 'MenuBar', 'None' );


%==============  GUI setup


MD.hFrameText = uicontrol('Parent', hFig, ...
    'Position',[ 5 5 70 14 ], ...
    'Style', 'text', 'FontName', 'FixedWidth', 'FontSize', 10, ...
    'String', ['1 / ', num2str(MD.nbFrames)] );

% setup frame Slider
MD.hFrame = uicontrol('Parent',hFig, ...
    'Position',[ 4 22 240 22 ], ...
    'Callback', [mfilename ' frameslider'], ...
    'Min',1, 'Max',MD.nbFrames, ...
    'SliderStep',[1/MD.nbFrames 5/MD.nbFrames], ...
    'Style','slider', ...
    'Value', 1 );

% setup STOP button
MD.hPlay = uicontrol('Parent',hFig, ...
    'Position',[ 80 2 50 20 ], ...
    'Callback', [mfilename ' play'], ...
    'FontWeight','Bold', ...
    'String','Play', ...
    'Style','pushbutton');

% setup Loop button
MD.hLoop = uicontrol('Parent',hFig, ...
    'Position',[ 135 2 50 20], ...
    'Callback', [mfilename ' loop'], ...
    'FontWeight','Bold', ...
    'String','Loop', ...
    'Style','pushbutton');

% setup Faster Button
MD.hFaster = uicontrol('Parent',hFig, ...
    'Position',[ 190 2 50 20], ...
    'Callback', [mfilename ' faster'], ...
    'FontWeight','Bold', ...
    'String','Faster', ...
    'Style','pushbutton');
    
% color min slider:
MD.hCmin = uicontrol('Parent',hFig, ...
    'Position',[ 290 2 128 20], ...
    'Callback', [mfilename ' Cmin'], ...
    'Min', 0, 'Max', 4096, ...
    'SliderStep',[1/2048 1/128], ...
    'Style','slider', ...
    'Value', MD.Cmin );

MD.hCminText = uicontrol('Parent', hFig, ...
    'Position',[ 250 4 40 16], ...
    'Style', 'text', 'FontName', 'FixedWidth', 'FontSize', 10, ...
    'ForegroundColor', 'b',...
    'String', num2str(MD.Cmin) );

% color max slider:
MD.hCmax = uicontrol('Parent',hFig, ...
    'Position',[ 290 24 128 20], ...
    'Callback', [mfilename ' Cmax'], ...
    'Min', 0, 'Max', 4096, ...
    'SliderStep',[1/2048 1/128], ...
    'Style','slider', ...
    'Value', MD.Cmax );

MD.hCmaxText = uicontrol('Parent', hFig, ...
    'Position',[ 250 26 40 16], ...
    'Style', 'text', 'FontName', 'FixedWidth', 'FontSize', 10, ...
    'ForegroundColor', 'b',...
    'String', num2str(MD.Cmax) );

%======================== show the first image:

axes('Units','pixels', 'Position',[0 50 figsize ]);

im = uint8( ( double( MD.movie(MD.frame).data ) - MD.Cmin ) * MD.Cscale );
MD.hImage = image( im, 'EraseMode', 'none' );
colormap( gray(256) );
set(gcf, 'DoubleBuffer','on');
axis off;
hold on;

%======================== handle mouse actions on the picture:

figid = gcf;
set( figid, 'WindowButtonDownFcn',   [mfilename ' mouse_down']);
set( figid, 'WindowButtonUpFcn',     [mfilename ' mouse_up']);
set( figid, 'WindowButtonMotionFcn', [mfilename ' mouse_motion']);
set( figid, 'KeyPressFcn',           [mfilename ' key_pressed']);
set( figid, 'DeleteFcn',             [mfilename ' stop']);

return;

%============================================================================================

%============================================================================================

function showImage(MD)

set( MD.hImage, 'CData', uint8( ( double( MD.movie(MD.frame).data ) - MD.Cmin ) * MD.Cscale ) );
set( MD.hFrame, 'Value', MD.frame );
set( MD.hFrameText, 'String', sprintf('%i / %i', MD.frame, MD.nbFrames));

return;