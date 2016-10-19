function [res, qu] = bindload(name);

% [res, qu] = bindload(name);

start = findstr(name, '\');
if isempty(start) start = 0; end
filepath = name(1:start(end));
start = length(name(start(end)+1:end));
res = struct('conc',0,'lam',0,'int',0);
lst = feval(@dir, [name '*']);
disp('Bla');
for j=1:length(lst) 
   tmp = lst(j).name(start+1:end-4)
   stm = findstr(tmp, '-');
   res(j).conc = str2num(strrep(tmp(1:stm-1),',','.'))/str2num(strrep(tmp(stm+1:end),',','.'));
   tmp = feval(@load, [filepath lst(j).name]); 
   res(j).lam = tmp(:,1); 
   res(j).int = tmp(:,2); 
end

qu=[]; 
for j=1:length(lst)
   qu = [qu; [res(j).conc sum(res(j).int)]];
end;
[tmp,ind]=sort(qu(:,1));
qu=qu(ind,:);
res = res(ind);

%p = simplex('bindfun',[10*max(qu(:,1)) 1 0.1], [0 0 0], [], [], [], qu(:,1), qu(:,2));
