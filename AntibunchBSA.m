close all
clear all

ww = 0.22; % estimated beam waist
load ConfoAnisoSat1
kk = 4;

if 0
    cd 'c:\Joerg\Doc\Fcs\BSA\Anastasia\'

    names = dir('2008-02*.mat');
    for j=4:4 %1:length(names)
        load(names(j).name);
        a(1) = []; b(1) = [];
        a = mean([a(1:2:2*floor(end/2)); a(2:2:2*floor(end/2))]);
        b = mean([b(1:2:2*floor(end/2))'; b(2:2:2*floor(end/2))'])';        

        para = [2.8954 12.4990 4.2467 1.7024 36.3598 0.3];
        close; para = Simplex('AntibunchDelay',para,[0 0 3 0 0 0],[inf inf 5 inf 100 inf],[],[],0.032*a,b,aniso(k),[],1); para
        for k=1:3
            close; para = Simplex('AntibunchDelay',para,para.*[1 1 1 0 0]',para.*[1 1 1 inf inf]',[],[],0.032*a,b,aniso,al,1); para
            close; para = Simplex('AntibunchDelay',para,[0 0 3 0 0],[inf inf 5 inf 100],[],[],0.032*a,b,aniso,al,1); para
        end
        %                res(kw,:) = para;
        %             end

        ax = axis;
        axis([-60 60 ax(3:4)])
        xlabel('time (ns)');
        ylabel('autocorrelation (a.u.)')

        text(20,ax(3) + 0.05*(ax(4)-ax(3)), ['\tau_{\itrot\rm} = ' mnum2str(para(5),2,1) ' ns'])

        %eval(['print -dpng -r300 ''' names(j).name(1:end-4) '_A''']);
    end
end

if 1
    %load('c:\Users\Public\Transfer\Antibunching\BSA_Alexa_antibunching_193C.mat')
    load c:\Users\Public\Transfer\Antibunching\antibunching_BSA_Atto647N_1pc_198Cbegin.mat

    autotime(1) = []; auto(1,:) = [];
    a = autoresolution*[-autotime(end:-1:1) autotime]';
    b = [[auto(end:-1:1,2); auto(1:end,1)], [auto(end:-1:1,4); auto(1:end,3)]];
    b = b(abs(a)<400,:); a = a(abs(a)<400);

    close; 
    para0 = [3.3 syncperiod 4 2 0.1];
    para0 = Simplex('AntibunchDelay',para0,[0 syncperiod 3 0 0],[inf syncperiod 5  inf inf],[],[],a(abs(a)<100),b(abs(a)<100,2),aniso(kk),1,1);
    [err,c,z] = AntibunchDelay(para0,a,b(:,2),aniso(kk),1,1);
    
    para = [para0(1:4)' 100];
    close;
    para=Simplex('AntibunchDelay0',para,[para0(1:3)' 0 0],[para0(1:3)' inf inf],[],[],a(abs(a)<200),b(abs(a)<200,:),z(abs(a)<200)-c(1),1); 
    para

    %     close; para = Simplex('AntibunchDelay',para,[0 syncperiod 3 0 0 0],[inf syncperiod 5 inf inf inf],[],[],a,b,aniso(kk),[],1); para
    %     kk = 10;
    %     for k=1:3
    %         close; para = Simplex('AntibunchDelay',para,[],[],[],[],a,b,aniso(kk),[],1); para
    %         close; para = Simplex('AntibunchDelay',para,[],[],[],[],a,b,aniso(kk),[],1); para
    %     end
    %
    %     ax = axis;
    %     axis([-60 60 ax(3:4)])
    %     xlabel('time (ns)');
    %     ylabel('autocorrelation (a.u.)')
    %
    %     text(20,ax(3) + 0.05*(ax(4)-ax(3)), ['\tau_{\itrot\rm} = ' mnum2str(para(5),2,1) ' ns'])

end

if 0
    cd 'c:\Joerg\Doc\Fcs\BSA\Anastasia\Translational\'

    names = dir('2008-0*.mat');
    for j=11:length(names)

        if(~isempty(findstr(names(j).name,'Atto655')))
            global pd; 
            pd = [0.1]; 
            res(j).dc = FCSFit(names(j).name,[500 150],0,1);
        else
            global pd; 
            pd = [0.1 10 20]; 
            res(j).dc = FCSFit(names(j).name,[500 150],2,1,[],[0 5 5; inf inf inf]);
        end
        
        eval(['print -dpng -r300 ''' names(j).name(1:end-4) '''']);
    end
end

return

loglog(1e3*w0,abs(pp./(max(po,[],2)*[1 1 1])),'r',1e3*w0,abs(po./(max(po,[],2)*[1 1 1])),'b',1e3*w0,op./(max(po,[],2)*[1 1 1]),'g')
set(gca,'xtick',2e3*[0.1 0.2 0.3])
axis([200 640 1e-2 1.1])
tmp = flipud(get(gca,'children'));
legend(tmp(1:3:end),{' || \times ||',' || \times \perp','\perp \times ||'},4)
j=2;
set(tmp(j),'linestyle','--')
set(tmp(3+j),'linestyle','--')
set(tmp(6+j),'linestyle','--')
j=3;
set(tmp(j),'linestyle',':')
set(tmp(3+j),'linestyle',':')
set(tmp(6+j),'linestyle',':')
for j=1:3
    if pp(1,j)<0
        set(tmp(j),'marker','o')
    end
    if po(1,j)<0
        set(tmp(3+j),'marker','o')
    end
end
xlabel('beam waist \itw\rm_0 (nm)')
ylabel('absolute coefficient values')
    