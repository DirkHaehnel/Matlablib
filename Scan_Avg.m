function [S] = Scan_Avg(name);

%  Returns sum of sscans of the same area. 

% (c) Ingo Gregor 2004

S = [0];

fname=strcat(name,'_01_*.t3r');
files = dir(fname);

for z=1:size(files,1)
	fname = files(1).name;
    zname = strrep(fname, name, '');
    zname = strrep(zname, '_01_', '_*_');
    zname = strcat(name,zname);
    zfiles = dir(zname);
    
    n_scans = size(zfiles,1);
    
	[tag, tim, head, tmp] = ScanRouterRead(zfiles(z).name);
	
%	A = tag(:,:,1)+tag(:,:,2);
	A = tag;
    
	if n_scans > 1 
       for r= 2:n_scans
           [tag, tim, head, tmp] = ScanRouterRead(zfiles(r).name);
%          B = tag(:,:,1)+tag(:,:,2);
           B = tag;
           A = ShiftMat(A,B);  
           mim(A, B, 'h'); drawnow; 
       end;
    end;
    if size(A,1)>size(S,1) 
        S(size(A,1), end,:) = 0;
    end;
    if size(A,2)>size(S,2) 
        S(end, size(A,2),:) = 0;
    end;
    
    S(:,:,z) = A;		

end;
    