% Pure longitudinal force MZ0
function [mz0] = MF96_MZ0(kappa, alpha, phi, Fz, tyre_data)

 % precode

  [alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(kappa, alpha, phi, Fz, tyre_data);

 % main code

  t = magic_formula_cos(alpha, alpha__t, Bt, Ct, Dt, Et);
  Mzr = magic_formula_cos(alpha, alpha__r, Br, 1, Dr, 0);
  mz0 = -t * MF96_FY0(kappa, alpha, phi, Fz, tyre_data) + Mzr;
  
 end
