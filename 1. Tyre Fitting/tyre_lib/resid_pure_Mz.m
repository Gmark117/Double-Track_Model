function res = resid_pure_Mz(fitting_params, P, MZ, KAPPA, ALPHA, GAMMA , FZ, tyre_data)
    % ----------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the Fx curve 
    %   Pacejka 1996 Magic Formula
    % ----------------------------------------------------------------------

    % Define MF coefficients
    tmp_tyre_data = tyre_data;
    for i = 1:length(fitting_params)
        tmp_tyre_data.(fitting_params{i}) = P(i);
    end 
    
    % Align Moment Equations
    res = 0;
    for i=1:length(ALPHA)
       mz0  = MF96_MZ0(KAPPA(i), ALPHA(i), GAMMA(i), FZ(i), tmp_tyre_data);
       res = res+(mz0-MZ(i))^2;
    end
    
    % Compute the residuals
    res = res/sum(MZ.^2);

end
