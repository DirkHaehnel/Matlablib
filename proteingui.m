% set of functions designed to implement a protein flight simulator
% Kathryn Lorenz, BioNB 441, Fall 2002

function proteingui(fcn)

if nargin == 0
   fcn = 'makeGUI';
end

switch fcn
    
%constructs "main" GUI 
case 'makeGUI'
    close all
    plotinfo.myname = mfilename;

    h0 = figure('FileName', 'protein_Gui.m', ...
    'MenuBar', 'none',...
    'Resize','off',...
    'NumberTitle','off',...
    'Name','Protein Navigator',...
    'Interruptible','off',...
    'Color',get(0,'DefaultUIControlBackgroundColor'));

    uicontrol('Parent', h0, ...
   'Units', 'points', ...
   'FontName', 'MS Sans Serif', ...
   'FontSize', 14, ...
   'FontWeight', 'bold', ...
   'ListboxTop', 0, ...
   'Position', [126.6 293 200 17.4], ...
   'String', 'Name of Protein', ...
   'Style', 'text', ...
   'Tag', 'TitleString');

plotinfo.ax = axes('Parent', h0, ...
   'Units', 'pixels', ...
   'CameraUpVector', [0 1 0], ...
   'Color', [1 1 1], ...
   'Position', [25 100 280 260], ...
   'Tag', 'Axes1', ...
   'XTick', 0, ...
   'YTick', 0, ...
   'ZTick', 0, ...
	'XColor', [0 0 0], ...
   'YColor', [0 0 0], ...
   'ZColor', [0 0 0]);

% Do axes contain a drawing?
plotinfo.proteinthere = 0;

% text labels

    uicontrol('Parent', h0, ...
   'Units', 'points', ...
   'Listboxtop',0, ...
   'Position', [285 245 70 17], ...
   'String', 'load protein:', ...
   'Style', 'text', ...
   'Tag', 'StaticText1');

    uicontrol('Parent', h0, ...
   'Units', 'points', ...
   'Listboxtop',0, ...
   'Position', [285 175 70 17], ...
   'String', 'display as:', ...
   'Style', 'text', ...
   'Tag', 'StaticText1');

    uicontrol('Parent', h0, ...
   'Units', 'points', ...
   'Listboxtop',0, ...
   'Position', [285 50 70 17], ...
   'String', 'Press ''Q'' to end', ...
   'Style', 'text', ...
   'Tag', 'StaticText1');

% browse button
    uicontrol('Parent', h0, ...
    'callback', [plotinfo.myname,' browse'], ... 
    'Units', 'points', ...
    'ListboxTop',0,...
    'Position', [285 230 70 17], ...
    'String', 'Browse', ...
    'Style', 'pushbutton', ...
    'Tag', 'Pushbutton1');

% start button

    uicontrol('Parent', h0, ...
    'Units', 'points', ...
    'callback', [plotinfo.myname,' start'], ...
    'ListboxTop',0,...
    'Position', [285 70 70 17], ...
    'String', 'Start', ...
    'Style', 'pushbutton', ...
    'Tag', 'Pushbutton1');

% radiobuttons

     uicontrol('Parent', h0, ...
    'Units', 'points', ...
    'Position', [285 120 70 55], ...
    'Style', 'frame', ...
    'Tag', 'Frame1');

    plotinfo.r1 = uicontrol('Parent', h0, ...
    'Units', 'points', ...
    'Position', [295 155 50 15], ...
    'String', 'line', ...
    'Style', 'radiobutton', ...
    'Tag', 'radiobutton1', ...
    'Value', 1,...
    'CallBack',[plotinfo.myname,' line'] );

    plotinfo.r2 = uicontrol('Parent', h0, ...
    'Units', 'points', ...
    'Position', [295 140 50 15], ...
    'String', 'ribbon', ...
    'Style', 'radiobutton', ...
    'Tag', 'radiobutton2', ...
    'Value', 0,...
    'CallBack',[plotinfo.myname,' ribbon'] );

    plotinfo.r3 = uicontrol('Parent', h0, ...
    'Units', 'points', ...
    'Position', [295 125 50 15], ...
    'String', 'ball & stick', ...
    'Style', 'radiobutton', ...
    'Tag', 'radiobutton3',...
    'Value', 0,...
    'CallBack',[plotinfo.myname,' ball'] );

   set(h0,'UserData',plotinfo);

case 'line'
    plotinfo=get(gcf,'UserData');  
    set(plotinfo.r1,'value',1);
    set(plotinfo.r2,'value',0);
    set(plotinfo.r3,'value',0);
    if plotinfo.proteinthere == 1
        cla
        plotinfo.ax = lineprot(plotinfo.file);
    end
    set(gcf,'userdata', plotinfo);
   
case 'ribbon'
    plotinfo=get(gcf,'UserData');  
    set(plotinfo.r1,'value',0);
    set(plotinfo.r2,'value',1);
    set(plotinfo.r3,'value',0);
    if plotinfo.proteinthere == 1
        cla
        plotinfo.ax = ribbon(plotinfo.file);
    end
    set(gcf,'userdata', plotinfo);
    
case 'ball'
    plotinfo=get(gcf,'UserData');  
    set(plotinfo.r1,'value',0);
    set(plotinfo.r2,'value',0);
    set(plotinfo.r3,'value',1);
    if plotinfo.proteinthere == 1
        cla
        plotinfo.ax=allAs(plotinfo.file);
    end
    set(gcf,'userdata', plotinfo);
    
case 'browse'
    cla
    plotinfo=get(gcf,'UserData'); 
    [filename pathname] = uigetfile('*.pdb', 'Open Protein File');
    plotinfo.file = filename;
    if get(plotinfo.r1,'value') ==1
        plotinfo.ax = lineprot(plotinfo.file);
    elseif get(plotinfo.r2, 'value') == 1
        plotinfo.ax = ribbon(plotinfo.file);
    else
        plotinfo.ax = allAs(plotinfo.file);
    end
    plotinfo.proteinthere =1;
    
    % a neat changing title for the main GUI
    newtitle = findobj(gcbf,'Tag','TitleString');
    set(newtitle,'String',plotinfo.file);
    set(gcf,'userdata', plotinfo);
    
% constructs "flight simulator" GUI    
case 'start'
    % these stored values come in handy again and again
    plotinfo=get(gcf,'UserData');
    plotinfo.plusangle=0;
    
    startfig = figure('FileName', 'flight.m', ...
    'MenuBar', 'none',...
    'Resize','off',...
    'NumberTitle','off',...
    'Name','Protein Flight Simulator',...
    'Interruptible','off',...
    'KeyPressFcn', 'simKeyEval(get(gcf,''CurrentCharacter''))',... % keystrokes
    'Color',[1 1 1],...
    'Position',[50 50 950 700]);


    startaxes = axes('Parent', startfig, ...
   'Units', 'pixels', ...
   'CameraPosition', [0 0 0], ...
   'CameraPositionMode', 'manual',...
   'CameraTarget', [20 20 0], ...
   'Color', [1 1 1], ...
   'Position', [2 2 949 699], ...
   'Tag', 'Axes1', ...
   'XTick', 0, ...
   'YTick', 0, ...
   'ZTick', 0, ...
	'XColor', [0 0 0], ...
   'YColor', [0 0 0], ...
   'ZColor', [0 0 0]);

    % for reset if necessary
    plotinfo.origcampos=get(gca,'cameraposition');
    plotinfo.origcamtar=get(gca,'cameratarget');
    plotinfo.origcamup=get(gca,'cameraupvector');
    
    if get(plotinfo.r1, 'value')==1
        startaxes=lineprot(plotinfo.file);
    elseif get(plotinfo.r2, 'value')==1
        startaxes=ribbon(plotinfo.file);
    else 
        startaxes=allAs(plotinfo.file);
            set(gca, 'cameratarget', [0.8219   22.4669   18.7550]);
            set(gca,'cameraposition',[ -159.3355 -186.2543  170.6485]);
            
    end
    set(gca, 'cameraviewangle', get(gca, 'cameraviewangle')/1.3),...
        'cameraupvector',[0 0 1];
     
    axis off
    hold on

    set(gcf,'UserData',plotinfo);

end
return


% plot line rendering of protein
function h11=lineprot(filename)

m=pickCA(filename);

splinex= spline(1:(size(m,1)),m(:,1),1:0.25:(size(m,1)));
spliney= spline(1:(size(m,1)),m(:,2),1:0.25:(size(m,1)));
splinez= spline(1:(size(m,1)),m(:,3),1:0.25:(size(m,1)));

for i = 1:(size(m,1)-1)
    verts{i}={[m(i,:);m(i+1,:)]};
end

h11 = plot3(splinex,spliney,splinez);


%plot ribbon rendering of protein
function l=ribbon(filename)

p=pickCA(filename);

splinex= spline(1:(size(p,1)),p(:,1),1:0.25:(size(p,1)));
spliney= spline(1:(size(p,1)),p(:,2),1:0.25:(size(p,1)));
splinez= spline(1:(size(p,1)),p(:,3),1:0.25:(size(p,1)));


verts = {[splinex' spliney' splinez']};

twistangle = {.1 * ones(size(verts{1},1),1)};
daspect([1 1 1])
l=streamribbon(verts,twistangle);
%-----Define viewing and lighting
axis tight

shading interp;
view(3);
camlight; lighting gouraud





% accepts keystrokes as commands during flight simulation
function simKeyEval(char)
plotinfo=get(gcf,'UserData');

pos=get(gca, 'cameraposition');
tar=get(gca,'cameratarget');
up=get(gca, 'cameraupvector')

switch char
    case 'q'
        quit_reply = questdlg('Really quit?');
        if(strcmp(quit_reply,'Yes'))
            close(gcf)
        end
    case 'r'        
        set(gca, 'cameraposition', plotinfo.origcampos);
        set(gca, 'cameratarget', plotinfo.origcamtar);
        set(gca, 'cameraupvector',plotinfo.origcamup);
    case num2str(5)
        vect1=tar-pos;
        vect1=3*vect1/norm(vect1);
        set(gca,'cameraposition',(pos+vect1));
    case num2str(0)
        vect1=tar-pos;
        vect1=vect1/norm(vect1);
        set(gca,'cameraposition',(pos-vect1));
    case num2str(2)
        vect1=tar-pos;
        nv=vect1*[cos(-pi/400) 0 -sin(-pi/400); 0 1 0; sin(-pi/400) 0 cos(-pi/400)];
        tar=nv+pos;
        set(gca,'cameratarget',tar);
    case num2str(8)
        vect1=tar-pos;
        nv=vect1*[cos(pi/400) 0 -sin(pi/400); 0 1 0; sin(pi/400) 0 cos(pi/400)];
        tar=nv+pos;
        set(gca,'cameratarget',tar);
    case num2str(4)
        vect1=tar-pos;
        nu=vect1*[cos(pi/400) sin(pi/400) 0; -sin(pi/400) cos(pi/400) 0; 0 0 1];
        tar=nu+pos;
        set(gca,'cameratarget',tar);
    case num2str(6)
        vect1=tar-pos;
        nu=vect1*[cos(-pi/400) sin(-pi/400) 0; -sin(-pi/400) cos(-pi/400) 0; 0 0 1];
        tar=nu+pos;
        set(gca,'cameratarget',tar);
    end
    
    
    function crd=allAs(pdbFileName)
crd=1; % error if it doesn't return something!

% comments aren't mine! but too many to erase
fid = fopen(pdbFileName);		% Open coordinate file
Eline = fgetl(fid);			% Skip till you get to ATOM lines
word = words(Eline); % word-break the line to words separated by space
% ~ is a logical not. Read next line if ~ TER and ~ ATOM
while (and(~strcmp(word{1},'TER'),~strcmp(word{1},'ATOM'))) 
	Fline = fgetl(fid);			% The first atom is N
	word = words(Fline);
end

i =0;
countN=1;
countC=1;
countO=1;
countS=1;
last2=[];

while ~strcmp(word{1},'TER');
    Fline = fgetl(fid);			% Find the CA of the next amino acid
	word = words(Fline);
    if or(strcmp(word{1},'TER'),strcmp(word{3},'OXT')) break; end;
    res=word{4};
       	%fprintf(fid_out,'%80s \n',Fline);
        switch  res
            case 'LYS'
                %N                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                %CA
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                %C
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                %O
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);           
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);

                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);               
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;             
                
            case 'VAL'
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;

                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);

                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);

                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last3=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last3;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
            case 'PHE'
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;  
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end                
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);

                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                              
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last3=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last3;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                verts={[last(1) last(2) last(3); last3(1) last3(2) last3(3)]};
                streamtube(verts,0.01,[1,10]);
            case 'GLY'
              
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end                
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
               
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);

                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                                                                
            case 'ARG'
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end                
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);  
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last3=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last3;
                countC=countC+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
            case 'CYS'

                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Svect{countS}=last;
                countS=countS+1;
            case 'GLU'

                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last3=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last3;
                countC=countC+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
            case 'LEU'

                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end                
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                 last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;              
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;            
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);              
                
                last3=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last3;
                countC=countC+1;             
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;             
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
            case 'ALA'
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;           
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;              
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;  
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);               
              
            case 'MET'

                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end                
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);               
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Svect{countS}=last;
                countS=countS+1;
                                             
                Fline = fgetl(fid);			 
		        word = words(Fline);
                                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);

                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
            case 'HIS'
               
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end                
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last3=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last3;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last4=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last4;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                verts={[last(1) last(2) last(3); last4(1) last4(2) last4(3)]};
                streamtube(verts,0.01,[1,10]);
            case 'TRP'
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end                
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);  
                %b
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
               Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                %g
                last3=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last3;
                countC=countC+1;
                
               Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last5=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last5;
                countC=countC+1;
                
               Fline = fgetl(fid);			 
		        word = words(Fline);
                
                last4=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last4;
                countC=countC+1;
                
               Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last4(1) last4(2) last4(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
               Nvect{countN}=last;
                countN=countN+1;
                
                verts={[last5(1) last5(2) last5(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
               Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);  
                
                verts={[last4(1) last4(2) last4(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                                last6=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last6;
                countC=countC+1;
                
               Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last6(1) last6(2) last6(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
               Fline = fgetl(fid);			 
		        word = words(Fline);                

                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
               Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last6(1) last6(2) last6(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                verts={[last4(1) last4(2) last4(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  

            case 'ASN'
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end                
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last3=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last3;
                countC=countC+1;
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
            case 'TYR'

                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end                
                
                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last3=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last3;
                countC=countC+1;
                
                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last4=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last4;
                countC=countC+1;
                
                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last4(1) last4(2) last4(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                                                Fline = fgetl(fid);			 
		        word = words(Fline);  
                
                verts={[last4(1) last4(2) last4(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                  
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
            case 'SER'

                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end                
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);

                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
            case 'THR'

                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end                
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);

                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);

                
                last3=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last3;
                countC=countC+1;
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
            case 'GLN'

                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end                
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);              
                
                last3=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last3;
                countC=countC+1;
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);

                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
            case 'ASP'

                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end                
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last3=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last3;
                countC=countC+1;
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
            case 'ILE'
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last;
                countO=countO+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                
                last3=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last3;
                countC=countC+1;
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last;
                countC=countC+1;
                
            case 'PRO'

                last=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Nvect{countN}=last;
                countN=countN+1;
                
                if length(last2)~=0
                    verts={[last2(1) last2(2) last2(3); last(1) last(2) last(3)]};
                    streamtube(verts,0.01,[1,10]);   
                end                
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last(1) last(2) last(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                ca=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=ca;
                countC=countC+1;
                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last2=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last2;
                countC=countC+1;
                                                
                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last2(1) last2(2) last2(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last4=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Ovect{countO}=last4;
                countO=countO+1;
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[ca(1) ca(2) ca(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);
                
                last3=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last3;
                countC=countC+1;
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline);
                
                verts={[last3(1) last3(2) last3(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);

                
                last5=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last5;
                countC=countC+1;
                
                                                                Fline = fgetl(fid);			 
		        word = words(Fline); 
                
                verts={[last5(1) last5(2) last5(3); str2double(word{6}) str2double(word{7}) str2double(word{8})]};
                streamtube(verts,0.01,[1,10]);                
                
                last5=[str2double(word{6}) str2double(word{7}) str2double(word{8})];
                Cvect{countC}=last5;
                countC=countC+1;
                
                verts={[last5(1) last5(2) last5(3); last(1) last(2) last(3)]};
                streamtube(verts,0.01,[1,10]);               
        end

    end
                for uu=1:length(Nvect)
                a=unitsphere(2);
                a=scale(a,0.5, 0.5, 0.5);
                a=translate(a,Nvect{uu}(1),Nvect{uu}(2),Nvect{uu}(3));
                a.facecolor='red';
                b=patch(a);
                set(b,'EdgeColor','none')
            end
            
            for u=1:length(Cvect)
                a=unitsphere(2);
                a=scale(a,0.5, 0.5, 0.5);
                a=translate(a,Cvect{u}(1),Cvect{u}(2),Cvect{u}(3));
                a.facecolor='green';
                b=patch(a);
                set(b,'EdgeColor','none')
            end
            
            for uuu = 1:length(Ovect)
                a=unitsphere(2);
                a=scale(a,0.5, 0.5, 0.5);
                a=translate(a,Ovect{uuu}(1),Ovect{uuu}(2),Ovect{uuu}(3));
                a.facecolor='blue';
                b=patch(a);
                set(b,'EdgeColor','none')
            end
            
            for uuuu=1:length(Svect)
                a=unitsphere(2);
                a=scale(a,0.5, 0.5, 0.5);
                a=translate(a,Svect{uuuu}(1),Svect{uuuu}(2),Svect{uuuu}(3));
                a.facecolor='green';
                b=patch(a);
                set(b,'EdgeColor','none')
            end


fclose(fid);

