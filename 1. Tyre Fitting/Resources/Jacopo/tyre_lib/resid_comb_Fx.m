function res = resid_comb_Fx(P,FX,KAPPA,ALPHA,GAMMA,FZ,tyre_data)

    % ------------------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the combined Fx curve 
    %  with Fz=Fz_nom, IA=0. Pacejka 1996 Magic Formula
    % ------------------------------------------------------------------------------

    % Define MF coefficients

    %Fz0 = 200*4.44822; % Nominal load 200 lbf
    
    tmp_tyre_data = tyre_data;
       
    tmp_tyre_data.rBx1 = P(1);
    tmp_tyre_data.rBx2 = P(2);
    tmp_tyre_data.rCx1 = P(3);
    tmp_tyre_data.rHx1 = P(4);
    
   %dfz = (Z - Fz0)./Fz0 ;
    
    % Longitudinal Force (Combined Longitudinal Slip) Equations
    res = 0;
    for i=1:length(KAPPA)
       fx  = MF96_FX(KAPPA(i), ALPHA(i), GAMMA, FZ, tmp_tyre_data);
       res = res+(fx-FX(i))^2;
       %res = res+(fx0/FX(i)-1)^2;
    end
    
    % Compute the residuals
    res = res/sum(FX.^2);
    

end



