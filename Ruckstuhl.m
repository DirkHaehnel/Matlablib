clear all
close all

if 0
    
    theta=(pi/10000:pi/5000:pi/2)';
    
    n0 = 1.5;
    n = n0 - 10.^(-(1:0.1:4));
    j = 1;
    [v,pc,ps] = DipoleL(theta,0,n0,n(j),n(j),[],0,[]);
    [v1,pc1,ps1] = DipoleL(theta,0,n(j),n(j),n0,[],0,[]);
    
    plot(sin(theta).*abs(v(:,1)).^2,-cos(theta).*abs(v(:,1)).^2,'r',...
        sin(theta).*sin(theta).^2*n0,-cos(theta).*sin(theta).^2*n0,'b:',...
        sin(theta).*abs(v1(:,1)).^2,cos(theta).*abs(v1(:,1)).^2,'r',...
        sin(theta).*sin(theta).^2*n0,cos(theta).*sin(theta).^2*n0,'b:',...
        -sin(theta).*abs(pc(:,1)).^2,-cos(theta).*abs(pc(:,1)).^2,'r',...
        -sin(theta).*cos(theta).^2*n0,-cos(theta).*cos(theta).^2*n0,'b:',...
        -sin(theta).*abs(pc1(:,1)).^2,cos(theta).*abs(pc1(:,1)).^2,'r',...
        -sin(theta).*cos(theta).^2*n0,cos(theta).*cos(theta).^2*n0,'b:',...
        -sin(theta).*abs(ps(:,1)).^2,-cos(theta).*abs(ps(:,1)).^2,'r',...
        -sin(theta).*n0,-cos(theta).*n0,'b:',...
        -sin(theta).*abs(ps1(:,1)).^2,cos(theta).*abs(ps1(:,1)).^2,'r',...
        -sin(theta)*n0,cos(theta)*n0,'b:');
    
    axis image
    ax = axis;
    ax = [-max(abs(ax(1:2))) max(abs(ax(1:2))) -max(abs(ax(3:4))) max(abs(ax(3:4)))];
    
    hold on
    plot(ax(1:2),[0 0],':g',[0 0],ax(3:4),':g');
    hold off
    axis(ax)
    
    text(0.6*ax(2),0.8*ax(4),['\itn\rm = ' num2str(n(j))])
    text(0.6*ax(2),-0.8*ax(4),['\itn\rm = ' num2str(n0)])    
    
    mov(j) = getframe;
    
    for j=2:length(n)
        [v,pc,ps] = DipoleL(theta,0,n0,n(j),n(j),[],0,[]);
        [v1,pc1,ps1] = DipoleL(theta,0,n(j),n(j),n0,[],0,[]);
        
        plot(sin(theta).*abs(v(:,1)).^2,-cos(theta).*abs(v(:,1)).^2,'r',...
            sin(theta).*sin(theta).^2*n0,-cos(theta).*sin(theta).^2*n0,'b:',...
            sin(theta).*abs(v1(:,1)).^2,cos(theta).*abs(v1(:,1)).^2,'r',...
            sin(theta).*sin(theta).^2*n0,cos(theta).*sin(theta).^2*n0,'b:',...
            -sin(theta).*abs(pc(:,1)).^2,-cos(theta).*abs(pc(:,1)).^2,'r',...
            -sin(theta).*cos(theta).^2*n0,-cos(theta).*cos(theta).^2*n0,'b:',...
            -sin(theta).*abs(pc1(:,1)).^2,cos(theta).*abs(pc1(:,1)).^2,'r',...
            -sin(theta).*cos(theta).^2*n0,cos(theta).*cos(theta).^2*n0,'b:',...
            -sin(theta).*abs(ps(:,1)).^2,-cos(theta).*abs(ps(:,1)).^2,'r',...
            -sin(theta).*n0,-cos(theta).*n0,'b:',...        
            -sin(theta).*abs(ps1(:,1)).^2,cos(theta).*abs(ps1(:,1)).^2,'r',...        
            -sin(theta)*n0,cos(theta)*n0,'b:',...        
            ax(1:2),[0 0],':g',[0 0],ax(3:4),':g');
        axis image
        axis(ax);
        
        text(0.6*ax(2),0.8*ax(4),['\itn\rm = ' num2str(n(j))])
        text(0.6*ax(2),-0.8*ax(4),['\itn\rm = ' num2str(n0)])    
        
        mov(j) = getframe;
    end
    
end

if 1 % Schwille paper

    ng = 1.5;
    nw = 1.33;
    zv = 2*pi*(0:0.25:1);
    
    % theta_crit = asin(nw/ng);
    theta=(pi/10000:pi/5000:pi/2)';

    for j=1:length(zv)
        [v(:,j),pc(:,j),ps(:,j)] = DipoleL(theta,zv(j),ng,nw,nw,[],zv(j),[]);
    end
    
    plot(theta/pi*180, (sin(theta)*ones(1,length(zv))).*(abs(v).^2+0.5*(abs(pc).^2+abs(ps).^2)));
    xlabel('angle (°)');
    ylabel('intensity (a.u.)')
    legend(num2str(zv'/2/pi),2)
end

