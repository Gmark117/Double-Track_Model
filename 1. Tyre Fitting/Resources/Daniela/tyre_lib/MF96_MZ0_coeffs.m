% Coefficients for Sel Aligning Torque 
function [alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(kappa, alpha, phi, Fz, tyre_data)

 % precode

  FZ0             = tyre_data.FZ0;
  pCy1            = tyre_data.pCy1;
  pDy1            = tyre_data.pDy1;
  pDy2            = tyre_data.pDy2;
  pDy3            = tyre_data.pDy3;
  pEy1            = tyre_data.pEy1;
  pEy2            = tyre_data.pEy2;
  pEy3            = tyre_data.pEy3;
  pEy4            = tyre_data.pEy4;
  pHy1            = tyre_data.pHy1;
  pHy2            = tyre_data.pHy2;
  pHy3            = tyre_data.pHy3;
  pKy1            = tyre_data.pKy1;
  pKy2            = tyre_data.pKy2;
  pKy3            = tyre_data.pKy3;
  pVy1            = tyre_data.pVy1;
  pVy2            = tyre_data.pVy2;
  pVy3            = tyre_data.pVy3;
  pVy4            = tyre_data.pVy4;

  R0              = tyre_data.R0;
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
  
  LCY             = tyre_data.LCY;
  LEY             = tyre_data.LEY;
  LFZ0            = tyre_data.LFZ0;
  LGAMMAY         = tyre_data.LGAMMAY;
  LHY             = tyre_data.LHY;
  LKA             = tyre_data.LKA;
  LKY             = tyre_data.LKY;
  LMR             = tyre_data.LMR;
  LMUY            = tyre_data.LMUY;
  LT              = tyre_data.LT;
  LVY             = tyre_data.LVY;
  

 % main code

  FZ01 = (LFZ0 * FZ0);
  dfz = Fz / FZ01 - 1;
  gamma__s = (phi * LGAMMAY);
  SHy = (dfz * pHy2 + pHy3 * gamma__s + pHy1) * LHY;
  SVy = Fz * (pVy1 + pVy2 * dfz + (dfz * pVy4 + pVy3) * gamma__s) * LVY * LMUY;
  alpha__y = alpha + SHy;
  Cy = pCy1 * LCY;
  mu__y = (dfz * pDy2 + pDy1) * (-pDy3 * gamma__s ^ 2 + 1) * LMUY;
  Dy = mu__y * Fz;
  Ey = (dfz * pEy2 + pEy1) * (1 - (pEy4 * gamma__s + pEy3) * Sign(alpha__y)) * LEY;
  Kya = FZ01 * pKy1 * sin(0.2e1 * atan((Fz / FZ01 / pKy2))) * (1 - pKy3 * abs(gamma__s)) * LFZ0 * LKA;
  By = Kya / Cy / Dy;
  gamma__z = (phi * LGAMMAY);
  SHf = SHy + SVy / Kya;
  SHt = qHz1 + qHz2 * dfz + (dfz * qHz4 + qHz3) * gamma__z;
  alpha__t = alpha + SHt;
  alpha__r = alpha + SHf;
  Bt = (dfz ^ 2 * qBz3 + dfz * qBz2 + qBz1) * (1 + qBz4 * gamma__z + qBz5 * abs(gamma__z)) * LKY / LMUY;
  Ct = qCz1;
  Dt = Fz * (dfz * qDz2 + qDz1) * (qDz4 * gamma__z ^ 2 + qDz3 * gamma__z + 1) * R0 / FZ0 * LT;
  Et = (dfz ^ 2 * qEz3 + dfz * qEz2 + qEz1) * (0.1e1 + (qEz5 * gamma__z + qEz4) * atan((Bt * Ct * alpha__t)));
  Br = (qBz9 * LKY / LMUY) + qBz10 * By * Cy;
  Dr = Fz * (qDz6 + qDz7 * dfz + (dfz * qDz9 + qDz8) * gamma__z) * R0 * LMUY * LMR;
  
 end
