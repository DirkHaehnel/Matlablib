function exc = DICExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over, focpos, pow, atf, resolution, maxm)

% important: DIC prism is always assumed to work along x-direction !!! 
% i.e. focus positions have to form a line along x!!!
% first focus is x-polarized, second focus is y-polarized

if isempty(d) || d<max(zfield(:))
    d = max(zfield(:));
end

if nargin<16 || isempty(resolution)
    resolution = 20;
end

if nargin<17 || isempty(maxm)
    maxm = 10;
end

if nargin<14 || isempty(pow)
    pow = 1;
end

if size(focpos,1)==1
    focpos = [focpos; focpos.*[1 -1 -1 1 1]];
end

if size(over,1)==1
    over = [1;1]*over;
end

ring = '';
[fxc1, fxs1, fyc1, fys1, fzc1, fzs1, rho, z] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over(1,:), focpos(1,:), atf, resolution, ring, maxm);
[fxc2, fxs2, fyc2, fys2, fzc2, fzs2] = GaussExc(rhofield, zfield, NA, fd, n0, n, n1, d0, d, d1, lamex, over(2,:), focpos(2,[1 3 2 4 5]).*[1 1 -1 1 1], atf, resolution, ring, maxm);
for j=1:maxm
    tmp = fxc2(:,:,j+1)*cos(j*pi/2) - fxs2(:,:,j)*sin(j*pi/2);
    fxs2(:,:,j) = fxc2(:,:,j+1)*sin(j*pi/2) + fxs2(:,:,j)*cos(j*pi/2);
    fxc2(:,:,j+1) = tmp;
    tmp = fyc2(:,:,j+1)*cos(j*pi/2) - fys2(:,:,j)*sin(j*pi/2);
    fys2(:,:,j) = fyc2(:,:,j+1)*sin(j*pi/2) + fys2(:,:,j)*cos(j*pi/2);
    fyc2(:,:,j+1) = tmp;
	tmp = -fyc2; fyc2 = fxc2; fxc2 = tmp; tmp = -fys2; fys2 = fxs2; fxs2 = tmp; 
    tmp = fzc2(:,:,j+1)*cos(j*pi/2) - fzs2(:,:,j)*sin(j*pi/2);
    fzs2(:,:,j) = fzc2(:,:,j+1)*sin(j*pi/2) + fzs2(:,:,j)*cos(j*pi/2);
    fzc2(:,:,j+1) = tmp;
end

exc.NA = NA;
exc.fd = fd;
exc.n0 = n0;
exc.n = n;
exc.n1 = n1;
exc.d0 = d0;
exc.d = d;
exc.d1 = d1;
exc.lambda = lamex;
exc.over = over;
exc.focpos = focpos;
exc.atf = atf;
exc.resolution = resolution;
exc.maxm = maxm;
exc.rho = rho;
exc.z = z;
exc.pow = pow;
exc.fxc1 = pow*fxc1;
exc.fxs1 = pow*fxs1;
exc.fyc1 = pow*fyc1;
exc.fys1 = pow*fys1;
exc.fzc1 = pow*fzc1;
exc.fzs1 = pow*fzs1;
exc.fxc2 = fxc2;
exc.fxs2 = fxs2;
exc.fyc2 = fyc2;
exc.fys2 = fys2;
exc.fzc2 = fzc2;
exc.fzs2 = fzs2;
