function S = Scan_Add(name);

%  Returns sum of scans of the same area. 

% (c) Ingo Gregor 2004

S = [0];

tmp = union(findstr(name,'\'),findstr(name,'/'));
path = name(1:tmp(end));
if ~isempty(tmp)
    name = name(tmp(end)+1:end);
end
files = dir([path name '_01_*.t3r']);

for z=1:size(files,1)
	fname = files(z).name
    zname = strrep(fname, name, '')
    zname = strrep(zname, '_01_', '_*_')
    zname = strcat(name,zname)
    zfiles = dir([path zname])
    
    n_scans = size(zfiles,1);
    
	[tag, tim, head, tmp] = ScanRouterRead([path zfiles(1).name]);
	zfiles(1).name
	A = sum(tag,3);
    
	if n_scans > 1 
       for r= 2:n_scans
           [tag, tim, head, tmp] = ScanRouterRead([path zfiles(r).name]);
           zfiles(r).name
           B = sum(tag,3);
           if size(A,1)>size(B,1) 
              B(size(A,1), end) = 0;
           end;
           if size(A,2)>size(B,2) 
              B(end, size(A,2)) = 0;
           end;
           if size(B,1)>size(A,1) 
              A(size(B,1), end) = 0;
           end;
           if size(B,2)>size(A,2) 
              A(end, size(B,2)) = 0;
           end;
           A = A + B;
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
    