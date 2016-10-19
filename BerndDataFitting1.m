if 1 % Rhodamin
    clear all
    
    dirs = {'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Rhodamin6G\2007-11-10\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Rhodamin6G\2007-11-11\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Rhodamin6G\2007-11-12\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Rhodamin6G\2007-11-12a\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Rhodamin6G\2007-11-13\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Rhodamin6G\2007-11-13a\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Rhodamin6G\2007-11-14\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Rhodamin6G\2007-11-14a\'};

    ind = {[2 3 4],...
        [1 4 6 7 8],...
        [1 2 3 4 5 6 7 8 9 10],...
        [3 5 6 8 9],...
        [1 2 3 4 5 7 10],...
        [3 5 6 7 8],...
        [3 4 5 6 8],...
        [3 5 6 7 9]};
    
    global pd
    for jd=1:1%:length(dirs)
        fnames = dir([dirs{jd} '*.fcsd']);
        clear data
        for jf=1:length(fnames)
            fin = fopen([dirs{jd} fnames(jf).name],'r');
            for j=1:196
                fgetl(fin);
            end
            y = fscanf(fin,'%f',inf);
            y = reshape(y,7,length(y)/7)';
            if ismember(jf,ind{jd})
                %                 semilogx(y(:,1),y(:,2:2:6))
                %                 pause
                %                 jf
                if ~exist('data')==1
                    data.y = y(3:end,2:2:6);
                    data.t = y(3:end,1)/1e3;
                else
                    data.y = data.y + y(3:end,2:2:6);
                end
            end
        end
        pd = [0.025 5];
        [dc v conc w0 a0 triplet c velo err] = FCSFit(data,[350 300],1,0,[100e3/60, [532 580]/1.33, 389],[[0 0];[inf 50]]);
        eval(['print -dpng -r300 ' dirs{jd} 'R6G_all.png'])
    end
end

if 0 % Atto465
    %clear all

    dirs = {'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Atto465\2007-11-07a\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Atto465\2007-11-19a\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Atto465\2007-11-20a\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Atto465\2007-11-21b\'};

    for jd=2:2; %1:length(dirs)
        fnames = dir([dirs{jd} '*.fcsd']);

        global pd
        pd = 1;
        clear data
        for jf=1:length(fnames)
            fin = fopen([dirs{jd} fnames(jf).name],'r');
            for j=1:196
                fgetl(fin);
            end
            y = fscanf(fin,'%f',inf);
            y = reshape(y,7,length(y)/7)';
            if ismember(jf,ind)
                if ~exist('data')==1
                    data.y = y(2:end,2:2:6);
                    data.t = y(2:end,1)/1e3;
                else
                    data.y = data.y + y(2:end,2:2:6);
                end
            end
            %data.y(data.t>1e-1,:) = [];
            %data.t(data.t>1e-1) = [];
            %[res(jd).dc(jf) res(jd).v(jf) res(jd).conc(jf) res(jd).w0(jf) res(jd).a0(jf) triplet c velo res(jd).err(jf)] = FCSFit(data,[300 100],0,0,[100e3/60, [470 510]/1.33, 370]);
        end
        [dc v conc w0 a0 triplet c velo err] = FCSFit(data,[300 100],0,0,[100e3/60, [470 510]/1.33, 370]);
        eval(['print -dpng -r300 ' dirs{jd} fnames(jf).name(1:end-5) '_b.png'])
    end

    %save Atto465 res
end


if 0 % Atto655
    clear all

    dirs = {'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Atto655\2007-10-26\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Atto655\2007-10-27\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Atto655\2007-10-28\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Atto655\2007-10-29\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Atto655\2007-10-29a\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Atto655\2007-10-30\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Atto655\2007-10-30b\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Atto655\2007-10-31\',...
        'c:\Joerg\Doc\Fcs\2Focus\Rhodamin\Atto655\2007-11-01\'};

    for jd=1:length(dirs)
        fnames = dir([dirs{jd} '*.fcsd']);

        global pd
        pd = 1;
        for jf=1:length(fnames)
            fin = fopen([dirs{jd} fnames(jf).name],'r');
            for j=1:196
                fgetl(fin);
            end
            y = fscanf(fin,'%f',inf);
            y = reshape(y,7,length(y)/7)';
            data.y = y(2:end,2:2:6);
            data.t = y(2:end,1)/1e3;
            %data.y(data.t>1e-1,:) = [];
            %data.t(data.t>1e-1) = [];

            [res(jd).dc(jf) res(jd).v(jf) res(jd).conc(jf) res(jd).w0(jf) res(jd).a0(jf) triplet res(jd).c(:,:,jf) velo res(jd).err(jf)] = FCSFit(data,[300 100],0,0,[100e3/60, [637 670]/1.33, 395.5]);
            eval(['print -dpng -r300 ' dirs{jd} fnames(jf).name(1:end-5) '.png'])
        end
    end

    save Atto655 res
end

if 0 % 1fFCS evaluation
    % rhodamine6g
    int = [1405 398 7429 3304 700 521 325 389]*0.0145;
    dc = [3.8 3.87 4.02 3.84 3.66 3.81 3.89 3.88];
    w0 = [334 334 345 339 315 329 330 336];
    [int,k] = sort(int); 
    dc = dc(k); w0 = w0(k);
    
    % atto655
    int = [20002 10182 5183 28035 53850 148274 1122 446 248]*0.0107;
    
    dc = {[5.06 5.12 5.07 5.01 5.00 5.02 5.04 4.95 5.02 5.06],...
        [4.66 4.7 4.96 5.03 4.83 4.92 4.9 4.66 4.91 4.8],...
        [4.74 4.57 4.73 4.71 4.75 4.66 4.71 4.2 4.56 4.48],...
        [4.92 5.14 5.25 5.09 5.09 5.15 5.11 4.78 4.9],...
        [5.13 5.09 5.06 5.05 4.85 5.02 5.05 4.86 5.08 4.88],...
        [5.17 5.2 5.11 5.12 5.1 5.05 4.9 5.03 5.07],...
        [4.15 4.21 4.17 4.27 4.16 4.35 4.07 4.21 4.19 4.16],...
        [4.49 4.12 3.98 4.17 4.75 4.16 4.35 4.29 3.99 4.14],...
        [4.41 4.59 4.27 3.94 3.78 4.47 4.07 4.71]};

    w0 = {[381 383 384 380 383 390 391 386 396 382],...
        [374 376 395 366 388 372 381 391 387 384],...
        [375 382 366 376 376 417 385 403 379 371],...
        [481 470 486 471 471 486 463 472 477],...
        [468 468 461 484 463 465 470 474 473 456],...
        [494 489 492 487 490 487 484 485 501],...
        [326 321 323 314 321 319 309 325 318 322],...
        [305 314 306 298 328 333 320 340 296 340],...
        [277 328 344 311 327 328 333 303]};
            
    %für Atto 465 (blau) ist der Umrechnungsfaktor int*0,0184 in uW
    int = [875 533 9516 1720]*0.0184;
    dc0 = [4.11 3.88 4.46 3.94];
    
    dc = {[4.31 4.58 4.61 4.25 4.71],...
        [4.36 3.96 4.24 4.28 4.4 4.15],...
        [4.4 4.47 4.45 4.49 4.52 4.4 4.29 4.47 4.42],...
        [4.48 4.28 4.69 4.25 4.4 4.47 4.46 4.28 4.61 4.36]};

    w0 = {[206 210 206 213 203],...
        [211 201 199 211 205 228],...
        [281 273 277 272 275 283 280 280 280],...
        [238 238 240 231 237 236 230 241 238 233]};
    
    
end
