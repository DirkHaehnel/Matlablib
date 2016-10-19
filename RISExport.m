% program RISExport

fin = fopen('D:\Joerg\SBData\Litera\Journals2.txt','r');
fout = fopen('D:\Joerg\SBData\Litera\Journals.txt','w');

tmpold = '';
while ~feof(fin)
    tmp = fgetl(fin);
    k = strfind(tmp,'Y2  - ');
    if ~isempty(k)
        kk = strfind(tmp,'.');
        if ~isempty(kk)
            yy = tmp(kk(2)+1:end);
            if str2num(yy)<8
                yy = ['20' yy];
            else
                yy = ['19' yy];
            end
            mm = tmp(kk(1)+1:kk(2)-1);
            if length(mm)<2
                mm = ['0' mm];
            end
            dd = tmp(7:kk(1)-1);
            if length(dd)<2
                dd = ['0' dd];
            end
            fprintf(fout,'%s\n', ['Y2  - ' yy '/' mm '/' dd '/']);
        else
            fprintf(fout,'%s\n', tmp);
        end
    else
        fprintf(fout,'%s\n', tmp);
    end
end

fclose(fin)
fclose(fout)

    