% Combined lateral/longitudinal forces FXY0
function [fxy0] = MF96_FXFYCOMB(kappa, alpha, phi, Fz, tyre_data)

 % precode

  [SHxa, Bxa, Cxa, Dxa, Gxa, SHyk, DVyk, SVyk, Byk, Cyk, Dyk, Gyk] = MF96_FXFYCOMB_coeffs(kappa, alpha, phi, Fz, tyre_data);

 % main code

  fxy0 = magic_formula(SHxa, Bxa, Cxa, Dxa, Gxa, SHyk, DVyk, SVyk, Byk, Cyk, Dyk, Gyk);
  
 end
