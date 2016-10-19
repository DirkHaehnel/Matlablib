function visc = Gdhcl2Visc(conc)

mol_ref =  [0 1.04 1.54 2.04 2.55 3.06 3.52 3.97 4.46 5.02 6.07]; 
visc_ref = [1.00253 1.03945 1.06426 1.09748 1.12893 1.16793 1.21812 1.27205 1.33496 1.43916 1.67286];

visc = interp1(mol_ref,visc_ref,conc,'cubic');


