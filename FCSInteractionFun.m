function [err, c, z, pd] = FCSInteractionFun(p,t,x,autotime,diffusion,rho,z,u)

diffusion = diffusion(:);
interact = FCSInteraction(p,autotime,rho,z,u);

pd = Simplex('Affinefit',10,0,[],[],[],t(:),x(:),autotime,[interact(:) diffusion], [], 1);

[err, c, z] = AffineFit(pd, t(:), x(:), autotime,[interact(:) diffusion], [], 1);

