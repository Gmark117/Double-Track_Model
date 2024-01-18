% Pure self alligning moment Mz0
% this function remap the scalar function to its vectorial form
function [mz0_vec] = MF96_MZ0_vec(kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

  mz0_vec = zeros(size(alpha_vec));
  for i = 1:length(alpha_vec)
   % precode
   [Br, Bt, Ct, Dr, Dt, Et, ~, ~, ~, ~, alpha__r, alpha__t, ~] = MF96_MZ0_coeffs(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
   mzr = Dr*cos(atan(Br*alpha__r))*cos(alpha_vec(i));
   t = Dt*cos(Ct*atan(Bt*alpha__t - Et*(Bt*alpha__t - atan(Bt*alpha__t))))*cos(alpha_vec(i));
   fy0 = MF96_FY0(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
   % main code 
   mz0_vec(i) = -fy0 * t + mzr;    % there's a possibility that it's  actually fy0 * t + mzr
  end
  
 end