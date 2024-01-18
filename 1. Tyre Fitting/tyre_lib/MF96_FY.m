% Combined slip longitudinal force FX
function [fy] = MF96_FY(kappa, alpha, phi, Fz, tyre_data)

 % precode

  [fy0] = MF96_FY0(kappa, alpha, phi, Fz, tyre_data); 
  [~, Gyk, SVyk] = MF96_FXcomb_coeffs(kappa, alpha, phi, Fz, tyre_data);

 % main code

  fy = Gyk * fy0 + SVyk;
  
 end
