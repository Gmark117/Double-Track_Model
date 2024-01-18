function res = resid_pure_Mz_varGamma(P,MZ,ALPHA,GAMMA,FZ,tyre_data)

    % ----------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the Mz curve 
    %  Pacejka 1996 Magic Formula
    % ----------------------------------------------------------------------

    % Define MF coefficients
    
    tmp_tyre_data = tyre_data;


    tmp_tyre_data.qBz4 = P(1);
    tmp_tyre_data.qBz5 = P(2);
    tmp_tyre_data.qDz8 = P(3);
    tmp_tyre_data.qDz9 = P(4);
    tmp_tyre_data.qEz5 = P(5);
    tmp_tyre_data.qHz3 = P(6);
    tmp_tyre_data.qHz4 = P(7);

    res = 0;
    for i=1:length(ALPHA)
       mz0  = MF96_MZ0(0, ALPHA(i), GAMMA(i), FZ, tmp_tyre_data);
       res = res+(mz0-MZ(i))^2;
    end
    
    % Compute the residuals
    res = res/sum(MZ.^2);

end

