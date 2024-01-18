% Implement the Magic formula cosine
function [mf_res] = magic_formula_cos(alpha, alpha_shifted, B, C, D, E)

 % precode

  

 % main code

  t1 = B * alpha_shifted;
  t2 = atan(t1);
  t6 = atan(-t1 + (t1 - t2) * E);
  t8 = cos(t6 * C);
  t10 = cos(alpha);
  mf_res = t10 * t8 * D;
  
 end
