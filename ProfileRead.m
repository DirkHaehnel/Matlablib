function [im,head] = ProfileRead(name);

if ~(name(end-2:end)=='wcb')
    disp('Not a wcb-file!');
   return
else
    fin = fopen(name,'r');
    if (fin==-1)
        errordlg('Cannot open specified file. Please try again.');
    else        
        head = struct('Signature', fread(fin, 1, 'uint32'));    
        head = setfield(head, 'Width', fread(fin, 1, 'uint32'));
        head = setfield(head, 'Height', fread(fin, 1, 'uint32'));
        head = setfield(head, 'BitsPerSample', fread(fin, 1, 'uint32'));
        head = setfield(head, 'Xpels', fread(fin, 1, 'uint32'));
        head = setfield(head, 'Ypels', fread(fin, 1, 'uint32'));
        
        fread(fin, 1, 'uint16');
        NumberOfPixels = head.Width*head.Height;       
        [y] = fread(fin, NumberOfPixels, 'uint16'); 
        im = reshape(y,head.Width,head.Height);
        fclose(fin);
    end
end

imagesc(im);
colorbar;