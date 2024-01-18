% Pure lateral force FY0
% this function remap the sclar function to its vectorial form
function [fy_vec] = MF96_FY_vec(kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

  
  fy_vec = zeros(size(alpha_vec));
  for i = 1:length(alpha_vec)
   fy_vec(i) = MF96_FY(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
  end
  
 end
