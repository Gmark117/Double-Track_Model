function res = resid_pure_Fy_varGamma(P,FY,ALPHA,GAMMA,FZ,tyre_data)

    % ----------------------------------------------------------------------
    %% Compute the residuals - least squares approach - to fit the Fy curve 
    %  with Fz=Fz_nom, IA=0. Pacejka 1996 Magic Formula
    % ----------------------------------------------------------------------

    % Define MF coefficients

    %Fz0 = 200*4.44822; % Nominal load 200 lbf
    
    tmp_tyre_data = tyre_data;
    % assigne computed parameter
%     tmp_tyre_data.pCy1 = P(1) ;
    tmp_tyre_data.pDy3 = P(1) ;
%     tmp_tyre_data.pDy2 = P(2) ;
%     tmp_tyre_data.pEy1 = P(4) ;
%     tmp_tyre_data.pEy2 = P(5) ;
    tmp_tyre_data.pEy3 = P(2) ;
    tmp_tyre_data.pEy4 = P(3) ;
%     tmp_tyre_data.pHy1 = P(8) ; 
%     tmp_tyre_data.pHy2 = P(9) ;
    tmp_tyre_data.pHy3 = P(4) ;
%     tmp_tyre_data.pKy1 = P(10) ;
%     tmp_tyre_data.pKy2 = P(11) ;
    tmp_tyre_data.pKy3 = P(5) ;
%     tmp_tyre_data.pVy1 = P(13) ;
%     tmp_tyre_data.pVy2 = P(14) ; 
    tmp_tyre_data.pVy3 = P(6) ;
    tmp_tyre_data.pVy4 = P(7) ;
    
    % Longitudinal Force (Pure Longitudinal Slip) Equations
    res = 0;
    for i=1:length(ALPHA)
       fy0  = MF96_FY0(ALPHA(i), 0, GAMMA(i), FZ, tmp_tyre_data);
       res = res+(fy0-FY(i))^2;
    end
    
    % Compute the residuals
    res = res/sum(FY.^2);

end

