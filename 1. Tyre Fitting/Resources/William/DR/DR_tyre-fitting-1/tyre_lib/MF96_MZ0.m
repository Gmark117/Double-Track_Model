% Pure longitudinal force MZ0
function [mz0] = MF96_MZ0(kappa, alpha, phi, Fz, tyre_data)

 % precode

  [gamma__z, SHf, SHt, alpha__t, alpha__r, Bt, Ct, Dt, Ed, Br, Dr] = MF96_MZ0_coeffs(kappa, alpha, phi, Fz, tyre_data);

 % main code

  mz0 = magic_formula(gamma__z, SHf, SHt, alpha__t, alpha__r, Bt, Ct, Dt, Ed, Br, Dr);
  
 end
