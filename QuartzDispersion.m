function [no, ne] = QuartzDispersion(lambda)

lambda = 1e-3*lambda;

Ao = 1.28604141;
Ae = 1.28851804;
Bo = 1.07044083;
Be = 1.09509924;
Co = 1.00585997e-2;
Ce = 1.02101864e-2;
Do = 1.10202242;
De = 1.15662475;
Fo = 100;
Fe = 100;

no = sqrt(Ao + Bo*lambda.^2./(lambda.^2-Co) + Do*lambda.^2./(lambda.^2-Fo));
ne = sqrt(Ae + Be*lambda.^2./(lambda.^2-Ce) + De*lambda.^2./(lambda.^2-Fe));