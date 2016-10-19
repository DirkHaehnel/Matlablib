% PhaseModAnalysis

if 1
    %     names = {'c:\Joerg\Doc\Fcs\2Focus\PhaseModulation\AT655 20uW Sin_10us.spc', ...
    %         'c:\Joerg\Doc\Fcs\2Focus\PhaseModulation\AT655 80uW Sin_10us again.spc', ...
    %         'c:\Joerg\Doc\Fcs\2Focus\PhaseModulation\AT655 80uW Sin_10us.spc', ...
    %         'c:\Joerg\Doc\Fcs\2Focus\PhaseModulation\AT655 200uW Sin_10us.spc', ...
    %         'c:\Joerg\Doc\Fcs\2Focus\PhaseModulation\AT655 160uW Sin_10us.spc'};

    %     names = {'c:\Joerg\Doc\Fcs\2Focus\PhaseModulation\AT655 10uW Sin_10us.spc', ...
    %         'c:\Joerg\Doc\Fcs\2Focus\PhaseModulation\AT655 100uW Sin_10us again.spc'};
    names = {'c:\Joerg\Doc\Fcs\2Focus\PhaseModulation\AT655 100uW Sin_10us.spc'};
    
    for j=1:length(names)
        name = names{j};
        [auto,autotime] = YouNew(name);

        eval(['save ''' name(1:end-4) ''' auto autotime'])

        ind = autotime>1e-6 & autotime<1e-1;
        autotime = autotime(ind);
        auto = auto(ind);

        close; 
        para = Simplex('PhaseModFcs',[430 130 10 4e2],[],[],[],[],100e3/60,[635 670]/1.33,433,1e6*autotime,auto);
        for k=1:2
            para = Simplex('PhaseModFcs',para,[],[],[],[],100e3/60,[635 670]/1.33,433,1e6*autotime,auto);
        end
        [err, c, z] = PhaseModFcs(para,100e3/60,[635 670]/1.33,433,1e6*autotime,auto);
        subplot(4,1,1:3); 
        semilogx(autotime,auto,autotime,z); 
        ylabel('correlation (a.u.)')
        subplot(4,1,4); 
        semilogx(autotime,auto-z,autotime,0*autotime,':')
        xlabel('time (s)'); ylabel('residuals')
        %ax = axis;
        %text(1e3, ax(3) + 0.8*(ax(4)-ax(3)), ['\itD\rm = ' mnum2str(para(4)*1e-2,1,2) '\cdot10^{-6} cm^2/s'])

        eval(['save ''' name(1:end-4) ''' para -append'])

        eval(['print -dpng -r300 ''' name(1:end-4) '_Full'''])
        
        if 0
            close
            ind = autotime>1e-1-1e-4;
            plot(autotime(ind),auto(ind),autotime(ind),z(ind));
            axis([1e-1-1e-4,1e-1,1850,2200])
            xlabel('time (s)'); ylabel('correlation (a.u.)')
        end
    end
end

if 0
    names = {'c:\Joerg\Doc\Fcs\2Focus\PhaseModulation\AT655 60uW Sin_10us.mat', ...
        'c:\Joerg\Doc\Fcs\2Focus\PhaseModulation\AT655 40uW Sin_10us.mat'};

    for j=1:length(names)
        name = names{j};
        load(name)

        ind = autotime>1e-6 & autotime<1e-1;
        autotime = autotime(ind);
        auto = auto(ind);

        close; para = Simplex('PhaseModFcs',[430 130 10 4e2],[],[],[],[],100e3/60,[635 670]/1.33,433,1e6*autotime,auto);
        xlabel('time (\mus)'); ylabel('correlation (a.u.)')
        ax = axis;
        text(1e3, ax(3) + 0.8*(ax(4)-ax(3)), ['\itD\rm = ' mnum2str(para(4)*1e-2,1,2) '\cdot10^{-6} cm^2/s'])

        %eval(['save ''' name ''' para -append'])

        eval(['print -dpng -r300 ''' name(1:end-4) '_Full'''])
    end
end

