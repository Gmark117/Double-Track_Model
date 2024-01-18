function res = resid_Fx_Comb(fitting_params, P, FX, KAPPA, ALPHA, GAMMA , FZ, tyre_data)
    % ----------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the Fx curve 
    %   Pacejka 1996 Magic Formula
    % ----------------------------------------------------------------------

    % Define MF coefficients
    tmp_tyre_data = tyre_data;
    for i = 1:length(fitting_params)
        tmp_tyre_data.(fitting_params{i}) = P(i);
    end 
    
    % Longitudinal Force (Pure Longitudinal Slip) Equations
    res = 0;
    for i=1:length(KAPPA)
       fx0  = MF96_FX(KAPPA(i), ALPHA(i), GAMMA(i), FZ(i), tmp_tyre_data);
       res = res+(fx0-FX(i))^2;
    end
    
    % Compute the residuals
    res = res/sum(FX.^2);

end

