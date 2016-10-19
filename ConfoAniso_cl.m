% Confocal anisotropy factors
close all
clear all

if 1 % compute with S0->S1 optical saturation 
    NA = 1.2;
    fd = 3e3;
    n0 = 1.333;
    n = 1.333;
    n1 = 1.333;
    d0 = [];
    d = [];
    d1 = [];
    lamex = 0.64;
    overv = 1e3:250:5e3;
    focpos = 30;
    av = 150/2; % confocal pinhole radius
    lamem = 0.670;
    mag = 60;
    zpin = 0;
    atf = [];
    resolution = 20;
    kappa = 1;

    rhofield = [0 10];
    zfield = [20 40];

    gammav = (15:15:90)/180*pi;
    % Kleine Vorschleife um alle notwendigen Parameter zu berechnen:
    for jo=1:length(overv)
        over{jo} = overv(jo);
        exc{jo} = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over{jo}, focpos, atf, resolution);
        ind{jo} = round(size(exc{jo}.z,2)/2);
        ff{jo} = abs(exc{jo}.fxc(:,ind{jo},1)).^2 + abs(exc{jo}.fyc(:,ind{jo},1)).^2 + abs(exc{jo}.fzc(:,ind{jo},1)).^2;
        tmp{jo} = Simplex('Gauss',[0 0.1],[0 0],[0 inf],[],[],[-flipud(exc{jo}.rho(:,1)); exc{jo}.rho(:,1)],[flipud(ff{jo}); ff{jo}]');
        w0(jo) = 2*tmp{jo}(2);
    end
    
    %Suchen nach dem jobmanager:
    jm = findResource('scheduler','type','jobmanager','Name','hal9000','LookupURL','skynet.physik3.gwdg.de')
    
    %Abhängige Files mitgeben:
    filedeps = cell(1,4);
    filedeps{1} = 'GaussExc2RotoFCS_cl.m';
    filedeps{2} = 'DipoleL.m';
    filedeps{3} = 'Fresnel.m';
    filedeps{4} = 'RotoDiffCoef.m';
    
    %Hier wird der Job gepackt
    job = createJob(jm);    
    set(job, 'FileDependencies', filedeps);
    
    for jo=1:length(overv)
        for jg = 1:length(gammav)
            % createTask(job, Funktionsname , Anz. Ausg.Param. , {Param. die überg. werden});
            createTask(job, @GaussExc2RotoFCS_cl, 5, {exc{jo}, lamem, mag, av, zpin, atf, gammav(jg)});

   %         [aniso(jo,jg).ux aniso(jo,jg).uy aniso(jo,jg).pp aniso(jo,jg).po aniso(jo,jg).op] = GaussExc2RotoFCS(exc{jo}, lamem, mag, av, zpin, atf, gamma);
        end
    end
    
    %Job abschicken:
    submit(job);
    tic;
    %Abwarten bis die Berechnung durchgelaufen ist:
    fprintf('Waiting for Jobdata\n');
    finished_jobs = 0;

    while finished_jobs ~= size(job, 2)
            if strcmp(get(job, 'State'), 'finished')
                [pending running finished] = findTask(job);
                without_error = 1;
                if without_error == 1 || size(finished, 2) == 0
                    finished_jobs = finished_jobs + 1;
                end
            end
        
        fprintf('%i Jobs of 1 finished\n', finished_jobs);
        if finished_jobs ~= size(job,2)
           pause(60) 
        end
    end
    
    %Jobs abholen:
        res = getAllOutputArguments(job);
        ind=1; 
        for jo=1:length(overv)
            for jg = 1:length(gammav)
                aniso(jo,jg).ux=res{ind,1}; 
                aniso(jo,jg).uy=res{ind,2}; 
                aniso(jo,jg).pp=res{ind,3}; 
                aniso(jo,jg).po=res{ind,4}; 
                aniso(jo,jg).op=res{ind,5}; 
                ind = ind+1; 
            end; 
        end
    
    %Speichern
        save ConfoAniso_Results_tst aniso w0 overv
 
    %Job vom Cluster löschen    
    
        destroy(job);
    
    %Ausrechnen wielange die Rechnung gedauert hat
    times = int32(toc);
    hh = num2str(fix(times/3600));
    min = num2str(fix((times-(fix(times/3600)*3600))/60));
    if size(min) == 1
       min = ['0' min]; 
    end
    sec = num2str(fix(times-((fix(times/3600)*3600)+(fix((times-(fix(times/3600)*3600))/60)*60))));
    if size(sec) == 1
       sec = ['0' sec]; 
    end 
    zeit = [num2str(hh) ':' num2str(min) ':' num2str(sec)];
    fprintf(['Calculation took: ' zeit ' hh:mm:ss' '\n']);        
end
