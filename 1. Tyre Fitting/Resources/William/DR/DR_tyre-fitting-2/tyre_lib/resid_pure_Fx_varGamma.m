function res = resid_pure_Fx_varGamma(P,FX,KAPPA,GAMMA,FZ,tyre_data)

    % ----------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the Fx curve 
    %  with Fz=Fz_nom, IA=0. Pacejka 1996 Magic Formula
    % ----------------------------------------------------------------------

    % Define MF coefficients

    %Fz0 = 200*4.44822; % Nominal load 200 lbf
    
    tmp_tyre_data = tyre_data;
    
    tmp_tyre_data.pDx3 = P(1);
    
    %dfz = (Z - Fz0)./Fz0 ;
    
    % Longitudinal Force (Pure Longitudinal Slip) Equations
    res = 0;
    for i=1:length(KAPPA)
       fx0  = MF96_FX0(KAPPA(i), 0, GAMMA(i), FZ, tmp_tyre_data);
       res = res+(fx0-FX(i))^2;
    end
    
    % Compute the residuals
    sumFX = sum(FX.^2);
    res = res/sumFX;

end
