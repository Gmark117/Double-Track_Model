% Coefficients for Magic Formula pure lateral force
function [gamma__z, SHf, SHt, alpha__t, alpha__r, Bt, Ct, Dt, Ed, Br, Dr] = MF96_MZ0_coeffs(kappa, alpha, phi, Fz, tyre_data)

 % precode

  qBz1            = tyre_data.qBz1;
  qBz10           = tyre_data.qBz10;
  qBz2            = tyre_data.qBz2;
  qBz3            = tyre_data.qBz3;
  qBz4            = tyre_data.qBz4;
  qBz5            = tyre_data.qBz5;
  qBz9            = tyre_data.qBz9;
  qCz1            = tyre_data.qCz1;
  qDz1            = tyre_data.qDz1;
  qDz2            = tyre_data.qDz2;
  qDz3            = tyre_data.qDz3;
  qDz4            = tyre_data.qDz4;
  qDz6            = tyre_data.qDz6;
  qDz7            = tyre_data.qDz7;
  qDz8            = tyre_data.qDz8;
  qDz9            = tyre_data.qDz9;
  qEz1            = tyre_data.qEz1;
  qEz2            = tyre_data.qEz2;
  qEz3            = tyre_data.qEz3;
  qEz4            = tyre_data.qEz4;
  qEz5            = tyre_data.qEz5;
  qHz1            = tyre_data.qHz1;
  qHz2            = tyre_data.qHz2;
  qHz3            = tyre_data.qHz3;
  qHz4            = tyre_data.qHz4;
  lambda__Fz0     = tyre_data.lambda__Fz0;
  lambda__Ky      = tyre_data.lambda__Ky;
  lambda__Mr      = tyre_data.lambda__Mr;
  lambda__mu__y   = tyre_data.lambda__mu__y;
  lambda__t       = tyre_data.lambda__t;
  lambda__gamma__y = tyre_data.lambda__gamma__y;
  

 % main code

  Fz01 = (lambda__Fz0 * Fz0);
  dfz = Fz / Fz01 - 1;
  gamma__z = (phi * lambda__gamma__y);
  SHf = (SHy + SVy / Kya);
  SHt = qHz1 + qHz2 * dfz + (dfz * qHz4 + qHz3) * gamma__z;
  alpha__t = alpha + SHt;
  alpha__r = alpha + SHf;
  Bt = (dfz ^ 2 * qBz3 + dfz * qBz2 + qBz1) * (1 + qBz4 * gamma__z + qBz5 * my_abs(gamma__z)) * lambda__Ky / lambda__mu__y;
  Ct = qCz1;
  Dt = Fz * (dfz * qDz2 + qDz1) * (qDz4 * gamma__z ^ 2 + qDz3 * gamma__z + 1) * R0 / Fz0 * lambda__t;
  Et = (dfz ^ 2 * qEz3 + dfz * qEz2 + qEz1) * (0.1e1 + (qEz5 * gamma__z + qEz4) * atan((Bt * Ct * alpha__t)));
  Br = qBz9 * lambda__Ky / lambda__mu__y + qBz10 * By * Cy;
  Dr = Fz * (qDz6 + qDz7 * dfz + (dfz * qDz9 + qDz8) * gamma__z) * R0 * lambda__mu__y * lambda__Mr;
  
 end
