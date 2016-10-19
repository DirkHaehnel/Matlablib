[buffer hbuffer] = tttrpro('f:\041111\Puffer_01.t3r','fimdafilter');
[gfp hgfp] = tttrpro('f:\041111\egfp_01.t3r','fimdafilter');
[p2x1 hp2x1] = tttrpro('f:\041111\P2X1_01.t3r','fimdafilter');
[glyr hglyr] = tttrpro('f:\041111\GlyR_01.t3r','fimdafilter');

nums = buffer.nums;
autotime = buffer.autotime;

fimda1 = cat(3,sum(buffer.fimda1(:,:,1:17),3), sum(gfp.fimda1(:,:,1:end-2),3), ...
        sum(glyr.fimda1(:,:,1:40),3), sum(p2x1.fimda1(:,:,1:18),3));

fimda2 = cat(3,sum(buffer.fimda2(:,:,1:17),3), sum(gfp.fimda2(:,:,1:end-2),3), ...
        sum(glyr.fimda2(:,:,1:40),3), sum(p2x1.fimda2(:,:,1:18),3));
