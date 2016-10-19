function [y, t] = FCSCrossReadRot(name,time)

if ischar(name)
    load(name);
else
    res = name;
end

dims = length(size(res.auto));

if nargin>1 && ~isempty(time)
    ind = res.autotime>=min(time) & res.autotime<=max(time);
else
    ind = res.autotime>0;
%     if (dims == 4)
%         ind = (res.autotime > 1e-8) & (res.autotime<3e-6);
%     end
end

if (dims == 6)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pulsed Rotational Diffusion Measurement

    t = res.autotime(ind);

    % 2 Laser -> 4 Detector

    % pol p : Detektor 1 + 3    Laser 2
    % pol o : Detektor 2 + 4    Laser 1

    % Curve 1:  pp x pp + oo x oo   
    y(:,1,:) = res.auto(ind,1,3,2,2,:) + res.auto(ind,3,1,2,2,:) + res.auto(ind,2,4,1,1,:) + res.auto(ind,4,2,1,1,:);

    % Curve 2:  pp x po + pp x op + po x pp + op x pp + ...
    % pppo
    pppo(:,1,:) = res.auto(ind,1,3,2,1,:) + res.auto(ind,3,1,2,1,:);
    % ppop
    pppo(:,2,:) = res.auto(ind,1,3,1,2,:) + res.auto(ind,3,1,1,2,:);
    % popp
    pppo(:,3,:) = res.auto(ind,1,2,2,2,:) + res.auto(ind,1,4,2,2,:) + res.auto(ind,3,2,2,2,:) + res.auto(ind,3,4,2,2,:);
    % oppp
    pppo(:,4,:) = res.auto(ind,2,1,2,2,:) + res.auto(ind,2,3,2,2,:) + res.auto(ind,4,1,2,2,:) + res.auto(ind,4,3,2,2,:);
    % ooop
    pppo(:,5,:) = res.auto(ind,2,4,1,2,:) + res.auto(ind,4,2,1,2,:);
    % oopo
    pppo(:,6,:) = res.auto(ind,2,4,2,1,:) + res.auto(ind,4,2,2,1,:);
    % opoo
    pppo(:,7,:) = res.auto(ind,2,1,1,1,:) + res.auto(ind,4,1,1,1,:) + res.auto(ind,2,3,1,1,:) + res.auto(ind,4,3,1,1,:);
    % pooo
    pppo(:,8,:) = res.auto(ind,1,2,1,1,:) + res.auto(ind,3,2,1,1,:) + res.auto(ind,1,4,1,1,:) + res.auto(ind,3,4,1,1,:);

    y(:,2,:) = sum(pppo, 2) / 5;

    % Curve 3:  po x po + op x op + po x op +...
    % ppoo
    ppoo(:,1,:) = res.auto(ind,1,3,1,1,:) + res.auto(ind,3,1,1,1,:);
    % oopp
    ppoo(:,2,:) = res.auto(ind,2,4,2,2,:) + res.auto(ind,4,2,2,2,:);
    % poop
    ppoo(:,3,:) = res.auto(ind,1,2,1,2,:) + res.auto(ind,1,4,1,2,:) + res.auto(ind,3,2,1,2,:) + res.auto(ind,3,4,1,2,:);
    % oppo
    ppoo(:,4,:) = res.auto(ind,2,1,2,1,:) + res.auto(ind,4,1,2,1,:) + res.auto(ind,2,3,2,1,:) + res.auto(ind,4,3,2,1,:);

    y(:,3,:) = sum(ppoo, 2) / 3;

    % Curve 4:  pp x oo + oo x pp 
    popo(:,1,:) = res.auto(ind,1,2,2,1,:) + res.auto(ind,2,1,1,2,:) + res.auto(ind,1,4,2,1,:) + res.auto(ind,4,1,1,2,:);
    popo(:,2,:) = res.auto(ind,3,2,2,1,:) + res.auto(ind,2,3,1,2,:) + res.auto(ind,3,4,2,1,:) + res.auto(ind,4,3,1,2,:);
    y(:,4,:) = sum(popo, 2) / 2;
elseif (dims == 4)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pulsed Rotational Diffusion Measurement 4 Dim
    
    t = res.autotime(ind)*1e9;
    
    L1 = 1;
    L2 = 2;
    D1 = 0;
    D2 = 2;
    D3 = 4;
    D4 = 6;    
    
    % pol p : Detektor 1 + 3    Laser 2
    % pol o : Detektor 2 + 4    Laser 1

    % Curve 1:  pp x pp + oo x oo   
    y(:,1,:) = res.auto(ind,D1+L2,D3+L2,:) + res.auto(ind,D3+L2,D1+L2,:) + res.auto(ind,D2+L1,D4+L1,:) + res.auto(ind,D4+L1,D2+L1,:);

    % Curve 2:  pp x po + pp x op + po x pp + op x pp + ...
    % pppo
    pppo(:,1,:) = res.auto(ind,D1+L2,D3+L1,:) + res.auto(ind,D3+L2,D1+L1,:);
    % ppop
    pppo(:,2,:) = res.auto(ind,D1+L1,D3+L2,:) + res.auto(ind,D3+L1,D1+L2,:);
    % popp
    pppo(:,3,:) = res.auto(ind,D1+L2,D2+L2,:) + res.auto(ind,D1+L2,D4+L2,:) + res.auto(ind,D3+L2,D2+L2,:) + res.auto(ind,D3+L2,D4+L2,:);
    % oppp
    pppo(:,4,:) = res.auto(ind,D2+L2,D1+L2,:) + res.auto(ind,D2+L2,D3+L2,:) + res.auto(ind,D4+L2,D1+L2,:) + res.auto(ind,D4+L2,D3+L2,:);
    % ooop
    pppo(:,5,:) = res.auto(ind,D2+L1,D4+L2,:) + res.auto(ind,D4+L1,D2+L2,:);
    % oopo
    pppo(:,6,:) = res.auto(ind,D2+L2,D4+L1,:) + res.auto(ind,D4+L2,D2+L1,:);
    % opoo
    pppo(:,7,:) = res.auto(ind,D2+L1,D1+L1,:) + res.auto(ind,D4+L1,D1+L1,:) + res.auto(ind,D2+L1,D3+L1,:) + res.auto(ind,D4+L1,D3+L1,:);
    % pooo
    pppo(:,8,:) = res.auto(ind,D1+L1,D2+L1,:) + res.auto(ind,D3+L1,D2+L1,:) + res.auto(ind,D1+L1,D4+L1,:) + res.auto(ind,D3+L1,D4+L1,:);

    y(:,2,:) = sum(pppo, 2) / 5;

    % Curve 3:  po x po + op x op + po x op +...
    % ppoo
    ppoo(:,1,:) = res.auto(ind,D1+L1,D3+L1,:) + res.auto(ind,D3+L1,D1+L1,:);
    % oopp
    ppoo(:,2,:) = res.auto(ind,D2+L2,D4+L2,:) + res.auto(ind,D4+L2,D2+L2,:);
    % poop
    ppoo(:,3,:) = res.auto(ind,D1+L1,D2+L2,:) + res.auto(ind,D1+L1,D4+L2,:) + res.auto(ind,D3+L1,D2+L2,:) + res.auto(ind,D3+L1,D4+L2,:);
    % oppo
    ppoo(:,4,:) = res.auto(ind,D2+L2,D1+L1,:) + res.auto(ind,D4+L2,D1+L1,:) + res.auto(ind,D2+L2,D3+L1,:) + res.auto(ind,D4+L2,D3+L1,:);

    y(:,3,:) = sum(ppoo, 2) / 3;

    % Curve 4:  pp x oo + oo x pp 
    popo(:,1,:) = res.auto(ind,D1+L2,D2+L1,:) + res.auto(ind,D2+L1,D1+L2,:) + res.auto(ind,D1+L2,D4+L1,:) + res.auto(ind,D4+L1,D1+L2,:);
    popo(:,2,:) = res.auto(ind,D3+L2,D2+L1,:) + res.auto(ind,D2+L1,D3+L2,:) + res.auto(ind,D3+L2,D4+L1,:) + res.auto(ind,D4+L1,D3+L2,:);
    y(:,4,:) = sum(popo, 2) / 2;    
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CW Rotational Diffusion Measurement

    t = res.autotime(ind) * 1e9;

    % 1 Laser CW -> 2 Detector

    % pol p : Detektor 1 + 3    Laser 2
    % pol o : Detektor 2 + 4    Laser 1

    % Curve 1:  pp x pp + oo x oo (cw: pp x pp only)
    y(:,1,:) = res.auto(ind,1,3,:) + res.auto(ind,3,1,:);

    % Curve 2:  pp x po + pp x op + po x pp + op x pp + ... (cw: op x pp + pp x op only)
    % popp
    pppo(:,1,:) = res.auto(ind,1,2,:) + res.auto(ind,1,4,:) + res.auto(ind,3,2,:) + res.auto(ind,3,4,:);
    % oppp
    pppo(:,2,:) = res.auto(ind,2,1,:) + res.auto(ind,2,3,:) + res.auto(ind,4,1,:) + res.auto(ind,4,3,:);
    
    y(:,2,:) = sum(pppo, 2) / 4;
    
    % Curve 3:  po x po + op x op + po x op +...  (cw: op x op only)
    y(:,3,:) = res.auto(ind,2,4,:) + res.auto(ind,4,2,:);
    
    % Curve 4:  pp x oo + oo x pp (cw: not possible)

end