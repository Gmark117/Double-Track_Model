function res = resid_Fy_Comb(fitting_params, P, FY, KAPPA, ALPHA, GAMMA , FZ, tyre_data)
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
    for i=1:length(ALPHA)
       fy0  = MF96_FY(KAPPA(i), ALPHA(i), GAMMA(i), FZ(i), tmp_tyre_data);
       res = res+(fy0-FY(i))^2;
    end
    
    % Compute the residuals
    res = res/sum(FY.^2);

end

