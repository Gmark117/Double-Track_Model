function res = resid_comb_Fy(P,FY,KAPPA,ALPHA,GAMMA,FZ,tyre_data)

    % ------------------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the combined Fy curve 
    %  Pacejka 1996 Magic Formula
    % ------------------------------------------------------------------------------

    % Define MF coefficients

    
    tmp_tyre_data = tyre_data;
       
    tmp_tyre_data.rBy1 = P(1);
    tmp_tyre_data.rBy2 = P(2);
    tmp_tyre_data.rBy3 = P(3);
    tmp_tyre_data.rCy1 = P(4);
    tmp_tyre_data.rHy1 = P(5);
    tmp_tyre_data.rVy1 = P(6);
    tmp_tyre_data.rVy4 = P(7);
    tmp_tyre_data.rVy5 = P(8);
    tmp_tyre_data.rVy6 = P(9);

    res = 0;
    for i=1:length(KAPPA)
       fy  = MF96_FY(KAPPA(i), ALPHA(i), GAMMA, FZ, tmp_tyre_data);
       res = res+(fy-FY(i))^2;
    end
    
    % Compute the residuals
    res = res/sum(FY.^2);
   
end