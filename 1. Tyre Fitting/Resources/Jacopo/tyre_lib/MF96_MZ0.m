% Pure self-aligning moment
function [mz0] = MF96_MZ0(kappa, alpha, phi, Fz, tyre_data)

 % precode

  [Br, Bt, Ct, Dr, Dt, Et, FZ01, SHf, SHt, dfz, alpha__r, alpha__t, gamma__z] = MF96_MZ0_coeffs(kappa, alpha, phi, Fz, tyre_data);
  mzr = Dr*cos(atan(Br*alpha__r))*cos(alpha);
  t = Dt*cos(Ct*atan(Bt*alpha__t - Et*(Bt*alpha__t - atan(Bt*alpha__t))))*cos(alpha);
  fy0 = MF96_FY0(kappa, alpha, phi, Fz, tyre_data);

 % main code

  mz0 = -fy0 * t + mzr;
  
 end
