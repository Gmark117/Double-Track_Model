function [Calfa_vec] = MF96_CorneringStiffness_lat(kappa_vec, alpha_vec, phi_vec, Fz_vec, tyre_data)

  
  Calfa_vec = zeros(size(alpha_vec));
  for i = 1:length(alpha_vec)
   % precode
   [alpha__y, By, Cy, Dy, Ey, SVy, ~, ~, ~] = MF96_FY0_coeffs(kappa_vec(i), alpha_vec(i), phi_vec(i), Fz_vec(i), tyre_data);
   % main code
    Calfa_vec(i) = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
  end
  
 end