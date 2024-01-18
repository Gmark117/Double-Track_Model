%% Initialisation
% 
% 
clc
clearvars 
close all   

% Set LaTeX as default interpreter for axis labels, ticks and legends
set(0,'defaulttextinterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultAxesFontSize',  16)
set(0,'DefaultLegendFontSize',16)

addpath('tyre_lib/')


to_rad = pi/180;
to_deg = 180/pi;

%% Select tyre dataset
%dataset path
data_set_path = 'dataset/';
% dataset selection and loading

data_set_LAT = 'Hoosier_B1464run23'; % pure lateral forces
data_set_LONG = 'Hoosier_B1464run30';  % braking/traction (pure log. force) + combined

% tyre geometric data:
% Hoosier	18.0x6.0-10
% 18 diameter in inches
% 6.0 section width in inches
% tread width in inches
diameter = 18*2.56; %
Fz0 = 700;   % [N] nominal load is given
R0  = diameter/2/100; % [m] get from nominal load R0 (m) *** TO BE CHANGED ***

fprintf('Loading dataset ...')
% LAT dataset
load ([data_set_path, 'Hoosier_B1464run23.mat']); % pure lateral
cut_start_LAT = 27760;
cut_end_LAT   = 54500;
% LONG dataset
data_LONG = load ([data_set_path, 'Hoosier_B1464run30.mat']); % pure longitudinal
cut_start_LONG = 19028;
cut_end_LONG   = 37643;

% select dataset portion
smpl_range_LAT = cut_start_LAT:cut_end_LAT;
smpl_range_LONG = cut_start_LONG:cut_end_LONG;

fprintf('completed!\n')

%% Plot raw data

figure(1)
tiledlayout(6,1)

ax_list(1) = nexttile; y_range = [min(min(-FZ),0) round(max(-FZ)*1.1)];
plot(-FZ)
hold on
plot([cut_start_LAT cut_start_LAT],y_range,'--r')
plot([cut_end_LAT cut_end_LAT],y_range,'--r')
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')

ax_list(2) = nexttile; y_range = [min(min(IA),0) round(max(IA)*1.1)];
plot(IA)
hold on
plot([cut_start_LAT cut_start_LAT],y_range,'--r')
plot([cut_end_LAT cut_end_LAT],y_range,'--r')
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(3) = nexttile; y_range = [min(min(SA),0) round(max(SA)*1.1)];
plot(SA)
hold on
plot([cut_start_LAT cut_start_LAT],y_range,'--r')
plot([cut_end_LAT cut_end_LAT],y_range,'--r')
title('Side slip')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(4) = nexttile; y_range = [min(min(SL),0) round(max(SL)*1.1)];
plot(SL)
hold on
plot([cut_start_LAT cut_start_LAT],y_range,'--r')
plot([cut_end_LAT cut_end_LAT],y_range,'--r')
title('Longitudinal slip')
xlabel('Samples [-]')
ylabel('[-]')

ax_list(5) = nexttile; y_range = [min(min(P),0) round(max(P)*1.1)];
plot(P)
hold on
plot([cut_start_LAT cut_start_LAT],y_range,'--r')
plot([cut_end_LAT cut_end_LAT],y_range,'--r')
title('Tyre pressure')
xlabel('Samples [-]')
ylabel('[psi]')

ax_list(6) = nexttile;  y_range = [min(min(TSTC),0) round(max(TSTC)*1.1)];
plot(TSTC,'DisplayName','Center')
hold on
plot(TSTI,'DisplayName','Internal')
plot(TSTO,'DisplayName','Outboard')
hold on
plot([cut_start_LAT cut_start_LAT],y_range,'--r')
plot([cut_end_LAT cut_end_LAT],y_range,'--r')
title('Tyre temperatures')
xlabel('Samples [-]')
ylabel('[degC]')

linkaxes(ax_list,'x')


%plot(SA,FY)

%% Select some specific data
% Cut crappy data and select only 12 psi data

vec_samples = 1:1:length(smpl_range_LAT);
tyre_data = table(); % create empty table

% store raw data in table
tyre_data.SL =  SL(smpl_range_LAT);
tyre_data.SA =  SA(smpl_range_LAT)*to_rad;
tyre_data.FZ =  -FZ(smpl_range_LAT);  % 0.453592  lb/kg
tyre_data.FX =  FX(smpl_range_LAT);
tyre_data.FY =  -FY(smpl_range_LAT);
tyre_data.MZ =  -MZ(smpl_range_LAT);
tyre_data.IA =  IA(smpl_range_LAT)*to_rad;

% Extract points at constant inclination angle
% To get points where the camber angle is almost 0, we use a tolerance
% Gamma_tol. MATLAB will look for data indexes that respect that condition.
% Then we populate Gamma_i with those indexes. This is done for different
% conditions and loads
GAMMA_tol = 0.05*to_rad;
idx.GAMMA_0 = 0.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 0.0*to_rad+GAMMA_tol;
idx.GAMMA_1 = 1.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 1.0*to_rad+GAMMA_tol;
idx.GAMMA_2 = 2.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 2.0*to_rad+GAMMA_tol;
idx.GAMMA_3 = 3.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 3.0*to_rad+GAMMA_tol;
idx.GAMMA_4 = 4.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 4.0*to_rad+GAMMA_tol;
idx.GAMMA_5 = 5.0*to_rad-GAMMA_tol < tyre_data.IA & tyre_data.IA < 5.0*to_rad+GAMMA_tol;
GAMMA_0  = tyre_data( idx.GAMMA_0, : );
GAMMA_1  = tyre_data( idx.GAMMA_1, : );
GAMMA_2  = tyre_data( idx.GAMMA_2, : );
GAMMA_3  = tyre_data( idx.GAMMA_3, : );
GAMMA_4  = tyre_data( idx.GAMMA_4, : );
GAMMA_5  = tyre_data( idx.GAMMA_5, : );

% Extract points at constant vertical load
% 200, 150, 50, 250, 100
FZ_tol = 100;
idx.FZ_220  = 220-FZ_tol < tyre_data.FZ & tyre_data.FZ < 220+FZ_tol;
idx.FZ_440  = 440-FZ_tol < tyre_data.FZ & tyre_data.FZ < 440+FZ_tol;
idx.FZ_700  = 700-FZ_tol < tyre_data.FZ & tyre_data.FZ < 700+FZ_tol;
idx.FZ_900  = 900-FZ_tol < tyre_data.FZ & tyre_data.FZ < 900+FZ_tol;
idx.FZ_1120 = 1120-FZ_tol < tyre_data.FZ & tyre_data.FZ < 1120+FZ_tol;
FZ_220  = tyre_data( idx.FZ_220, : );
FZ_440  = tyre_data( idx.FZ_440, : );
FZ_700  = tyre_data( idx.FZ_700, : );
FZ_900  = tyre_data( idx.FZ_900, : );
FZ_1120 = tyre_data( idx.FZ_1120, : );

% The slip angle is varied continuously between -4 and +12° and then
% between -12° and +4° for the pure slip case

% The slip angle is varied step wise for longitudinal slip tests
% 0° , - 3° , -6 °
% SA_tol = 0.5*to_rad;
% idx.SA_0    =  0-SA_tol            < tyre_data.SA & tyre_data.SA < 0+SA_tol;
% idx.SA_3pos =  (3*to_rad-SA_tol)   < tyre_data.SA & tyre_data.SA < 3*to_rad+SA_tol;
% idx.SA_4neg =  -(4*to_rad+SA_tol)  < tyre_data.SA & tyre_data.SA < -4*to_rad+SA_tol;
% idx.SA_12pos = (12*to_rad-SA_tol)  < tyre_data.SA & tyre_data.SA < 12*to_rad+SA_tol;
% idx.SA_12neg = -(12*to_rad+SA_tol) < tyre_data.SA & tyre_data.SA < -12*to_rad+SA_tol;
% SA_0     = tyre_data( idx.SA_0, : );
% SA_3pos  = tyre_data( idx.SA_3pos, : );
% SA_4neg  = tyre_data( idx.SA_4neg, : );
% SA_12pos  = tyre_data( idx.SA_12pos, : );
% SA_12neg  = tyre_data( idx.SA_12neg, : );

figure(2)
tiledlayout(3,1)

ax_list(1) = nexttile;
plot(tyre_data.IA*to_deg)
hold on
plot(vec_samples(idx.GAMMA_0),GAMMA_0.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_1),GAMMA_1.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_2),GAMMA_2.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_3),GAMMA_3.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_4),GAMMA_4.IA*to_deg,'.');
plot(vec_samples(idx.GAMMA_5),GAMMA_5.IA*to_deg,'.');
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(2) = nexttile;
plot(tyre_data.FZ)
hold on
plot(vec_samples(idx.FZ_220),FZ_220.FZ,'.');
plot(vec_samples(idx.FZ_440),FZ_440.FZ,'.');
plot(vec_samples(idx.FZ_700),FZ_700.FZ,'.');
plot(vec_samples(idx.FZ_900),FZ_900.FZ,'.');
plot(vec_samples(idx.FZ_1120),FZ_1120.FZ,'.');
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')


ax_list(3) = nexttile;
plot(tyre_data.SA*to_deg)
% hold on
% plot(vec_samples(idx.SA_0),   SA_0.SA*to_deg,'.');
% plot(vec_samples(idx.SA_3pos),SA_3pos.SA*to_deg,'.');
% plot(vec_samples(idx.SA_4neg),SA_4neg.SA*to_deg,'.');
% plot(vec_samples(idx.SA_12pos),SA_12pos.SA*to_deg,'.');
% plot(vec_samples(idx.SA_12neg),SA_12neg.SA*to_deg,'.');
title('Side slip')
xlabel('Samples [-]')
ylabel('[rad]')

%% Intersect tables to obtain specific sub-datasets

[TData0, ~] = intersect_table_data( tyre_data, GAMMA_0, FZ_700 );

%[tableSelectedData, ~] = intersectTableData( KAPPA_00, GAMMA_00, VX_10 );
%[TDataSub, ~] = intersectTableData( ALPHA_00, GAMMA_00, VX_10, FZ_6500 );

%[tableSelectedData, ~] = intersectTableData( KAPPA_00, GAMMA_00, VX_10, FZ_6500 );

% get data for tyre deflection (radius) versus speed
%[TDataSubRho, ~] = intersectTableData( KAPPA_00, ALPHA_00, GAMMA_00, FZ_6500 );

%% plot_selected_data

figure('Name','3. Selected-data')
plot_selected_data(TData0);

%% FITTING 
% Initialise tyre data
tyre_coeffs = initialise_tyre_data(R0, Fz0);

%% Fitting with Fz=Fz_nom= 700N and camber=0  alpha = 0 VX= 10
% ------------------
% lat slip

% Fit the coeffs {pCx1, pDx1, pEx1, pEx4, pKx1, pHx1, pVx1}
FZ0 = mean(TData0.FZ);

zeros_vec = zeros(size(TData0.SA));
ones_vec  = ones(size(TData0.SA));

FY0_guess = MF96_FY0_vec(zeros_vec, TData0.SA , zeros_vec, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);

% check guess 
figure(4)
plot(TData0.SA,TData0.FY,'.')
hold on
plot(TData0.SA,FY0_guess,'-')

% Plot raw data and initial guess
% figure()
% plot(TDataSub.KAPPA,TDataSub.FX,'o')
% hold on
% plot(TDataSub.KAPPA,FX0_guess,'x')

% Guess values for parameters to be optimised
%    [pCy1 pDy1 pEy1 pHy1 pKy1 pKy2 pVy1]
P0 = [  2,   2,   1,  0.1, 40,  70,  0.1];
% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
%    [pCy1 pDy1 pEy1 pHy1 pKy1 pKy2 pVy1]
lb = [-10,  -7,  -1,  -1,  -1,   0,   -1];
ub = [  0,   0,   1,   1,  50,  10,    1];

ALPHA_vec = TData0.SA;
FY_vec    = TData0.FY;

% check guess
SA_vec = -0.3:0.001:0.3;

% 
% figure
% plot(KAPPA_vec,FX_vec,'.')
% hold on
% plot(SL_vec,FX0_fz_nom_vec,'.')
% 


% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 

[P_fz_nom, ~] = fmincon(@(P)resid_pure_Fy(P,FY_vec, ALPHA_vec,0,FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre data with new optimal values                             
tyre_coeffs.pCy1 = P_fz_nom(1) ; % 1
tyre_coeffs.pDy1 = P_fz_nom(2) ;
% tyre_coeffs.pDy2 = P_fz_nom(3) ;
tyre_coeffs.pEy1 = P_fz_nom(3) ;
% tyre_coeffs.pEy2 = P_fz_nom(5) ;
% tyre_coeffs.pEy3 = P_fz_nom(6) ;
% tyre_coeffs.pEy4 = P_fz_nom(7) ;
tyre_coeffs.pHy1 = P_fz_nom(4) ; 
% tyre_coeffs.pHy2 = P_fz_nom(9) ;
tyre_coeffs.pKy1 = P_fz_nom(5) ;
tyre_coeffs.pKy2 = P_fz_nom(6) ;
% tyre_coeffs.pKy3 = P_fz_nom(12) ;
tyre_coeffs.pVy1 = P_fz_nom(7) ;
% tyre_coeffs.pVy2 = P_fz_nom(14) ;

FY0_fz_nom_vec = MF96_FY0_vec(zeros(size(SA_vec)) , SA_vec , zeros(size(SA_vec)), ...
                              FZ0.*ones(size(SA_vec)),tyre_coeffs);

figure('Name','5. Fy0(Fz0)')
plot(TData0.SA,TData0.FY,'o')
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
plot(SA_vec,FY0_fz_nom_vec,'-','LineWidth',2)
xlabel('$\alpha$ [-]')
ylabel('$F_{Y0}$ [N]')

res_Fy0=resid_pure_Fy(P_fz_nom,FY_vec, ALPHA_vec,0,tyre_coeffs.FZ0, tyre_coeffs);

fprintf('R-squared_FzNom = %6.3f\n',1-res_Fy0);

%% Fit coefficient with variable load
% extract data with variable load
[TDataDFz, ~] = intersect_table_data( tyre_data, GAMMA_0 );

% long slip

% Fit the coeffs {pCx1, pDx1, pEx1, pEx4, pKx1, pHx1, pVx1}
FZ0 = mean(TData0.FZ);

zeros_vec = zeros(size(TDataDFz.SA));
ones_vec  = ones(size(TDataDFz.SA));

% figure
% subplot(2,1,1)
% plot(TDataDFz.SA)
% subplot(2,1,2)
% plot(TDataDFz.FY)
% 
FY0_guess = MF96_FY0_vec(zeros_vec, TDataDFz.SA, zeros_vec, FZ0*ones_vec, tyre_coeffs);
% FY0_fz_nom_vec = MF96_FY0_vec(zeros(size(TDataDFz.SA)) , TDataDFz.SA , zeros(size(TDataDFz.SA)), ...
%                           667*ones(size(TDataDFz.SA)),tyre_coeffs);
% FY0_dfz_vec = MF96_FY0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), ...
%                            TDataDFz.FZ,tyre_coeffs);

% check guess 
figure(6)
plot(TDataDFz.SA,TDataDFz.FY,'.')
hold on
plot(TDataDFz.SA,FY0_guess,'-')
% plot(TDataDFz.SA,FY0_fz_nom_vec,'-','LineWidth',2)

% Plot raw data and initial guess
% figure()
% plot(TDataSub.ALPHA,TDataSub.FY,'o')
% hold on
% plot(TDataSub.ALPHA,FY0_guess,'y')

% Guess values for parameters to be optimised
%    [pCy1 pDy1 pDy2 pEy1 pEy2 pEy3 pEy4 pHy1 pHy2 pKy1 pKy2 pKy3 pVy1 pVy2]
P0 = [  1,   2,   1,  0,   0,   1,   0,   1,    0,   0,   1,   2,   0,   1]; 


% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%    [                                                                      |
lb = [-10, -7, -10, -20, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10];
ub = [  0,  0,  10,  35,  10,  35,  10,  10,  10,  35,  10,  10,  10,  10];


ALPHA_vec = TDataDFz.SA;
FY_vec    = TDataDFz.FY;
FZ_vec    = TDataDFz.FZ;

% check guess
SA_vec = -0.3:0.001:0.3;
FY0_dfz_vec = MF96_FY0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), ...
                           TDataDFz.FZ,tyre_coeffs);
% 
% figure
% plot(KAPPA_vec,FX_vec,'.')
% hold on
% plot(SL_vec,FX0_dfz_vec,'.')


% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 

[P_dfz,~,exitflag] = fmincon(@(P)resid_pure_Fy_varFz(P,FY_vec, ALPHA_vec,0,FZ_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

disp(exitflag)
% Change tyre data with new optimal values                             
tyre_coeffs.pCy1 = P_dfz(1) ; % 1
tyre_coeffs.pDy1 = P_dfz(2) ;
tyre_coeffs.pDy2 = P_dfz(3) ;
tyre_coeffs.pEy1 = P_dfz(4) ;
tyre_coeffs.pEy2 = P_dfz(5) ;
tyre_coeffs.pEy3 = P_dfz(6) ;
tyre_coeffs.pEy4 = P_dfz(7) ;
tyre_coeffs.pHy1 = P_dfz(8) ; 
tyre_coeffs.pHy2 = P_dfz(9) ;
tyre_coeffs.pKy1 = P_dfz(10) ;
tyre_coeffs.pKy2 = P_dfz(11) ;
tyre_coeffs.pKy3 = P_dfz(12) ;
tyre_coeffs.pVy1 = P_dfz(13) ;
tyre_coeffs.pVy2 = P_dfz(14) ;


res_FY0_dfz_vec = resid_pure_Fy_varFz(P_dfz,FY_vec,SA_vec,0 , FZ_vec,tyre_coeffs);

tmp_zeros = zeros(size(SA_vec));
tmp_ones = ones(size(SA_vec));


FY0_fz_var_vec1 = MF96_FY0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec2 = MF96_FY0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_440.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec3 = MF96_FY0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec4 = MF96_FY0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
FY0_fz_var_vec5 = MF96_FY0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);


figure('Name','7. Fy0(Fz0)')
plot(TDataDFz.SA,TDataDFz.FY,'o')
hold on
%plot(TDataSub.ALPHA,FY0_fz_nom_vec,'-')
%plot(SA_vec,FY0_dfz_vec,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec1,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec2,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec3,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec4,'-','LineWidth',2)
plot(SA_vec,FY0_fz_var_vec5,'-','LineWidth',2)

xlabel('$\alpha$ [-]')
ylabel('$F_{y0}$ [N]')


[alpha__y, By, Cy, Dy, Ey, SVy] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_220.FZ), tyre_coeffs);
Calfa_vec1_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_440.FZ), tyre_coeffs);
Calfa_vec2_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_700.FZ), tyre_coeffs);
Calfa_vec3_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_900.FZ), tyre_coeffs);
Calfa_vec4_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);
[alpha__y, By, Cy, Dy, Ey, SVy] =MF96_FY0_coeffs(0, 0, 0, mean(FZ_1120.FZ), tyre_coeffs);
Calfa_vec5_0 = magic_formula_stiffness(alpha__y, By, Cy, Dy, Ey, SVy);

Calfa_vec1 = MF96_CorneringStiffness_LAT(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec2 = MF96_CorneringStiffness_LAT(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_440.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec3 = MF96_CorneringStiffness_LAT(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec4 = MF96_CorneringStiffness_LAT(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
Calfa_vec5 = MF96_CorneringStiffness_LAT(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);

figure('Name','8. C_alpha')
subplot(2,1,1)
hold on
%plot(TDataSub.ALPHA,FY0_fz_nom_vec,'-')
plot(mean(FZ_220.FZ),Calfa_vec1_0,'+','LineWidth',2)
plot(mean(FZ_440.FZ),Calfa_vec2_0,'+','LineWidth',2)
plot(mean(FZ_700.FZ),Calfa_vec3_0,'+','LineWidth',2)
plot(mean(FZ_900.FZ),Calfa_vec4_0,'+','LineWidth',2)
plot(mean(FZ_1120.FZ),Calfa_vec5_0,'+','LineWidth',2)
legend({'$Fz_{220}$','$Fz_{440}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'})

subplot(2,1,2)
hold on
%plot(TDataSub.ALPHA,FY0_fz_nom_vec,'-')
plot(SA_vec,Calfa_vec1,'-','LineWidth',2)
plot(SA_vec,Calfa_vec2,'-','LineWidth',2)
plot(SA_vec,Calfa_vec3,'-','LineWidth',2)
plot(SA_vec,Calfa_vec4,'-','LineWidth',2)
plot(SA_vec,Calfa_vec5,'-','LineWidth',2)
legend({'$Fz_{220}$','$Fz_{440}$','$Fz_{700}$','$Fz_{900}$','$Fz_{1120}$'})

%% Fit coefficient with variable camber

% extract data with variable load
[TDataGamma, ~] = intersect_table_data( tyre_data, FZ_700 );

% Fit the coeffs { pDx3}

% Guess values for parameters to be optimised
P0 = [  1,   2,   1,  0,   0,   1,   0,   1,    0,   0,   1,   2,   0,   1]; 

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%lb = [0, 0,  0, 0,  0,  0,  0];
%ub = [2, 1e6,1, 1,1e1,1e2,1e2];
lb = [-10, -7, -10, -100, -10, -100, -100, -100, -100, -100, -100, -100, -100, -100];
ub = [  0,  0,  10,  100,  10,  100,  100,  100,  100,  100,  100,  100,  100,  100];


zeros_vec = zeros(size(TDataGamma.SA));
ones_vec  = ones(size(TDataGamma.SA));

ALPHA_vec = TDataGamma.SA;
GAMMA_vec = TDataGamma.IA; 
FY_vec    = TDataGamma.FY;
FZ_vec    = TDataGamma.FZ;

figure(9)
plot(ALPHA_vec,FY_vec);


% LSM_pure_Fy returns the residual, so minimize the residual varying Y. It
% is an unconstrained minimization problem 
[P_varGamma,fval,exitflag] = fmincon(@(P)resid_pure_Fy_varGamma(P,FY_vec, ALPHA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values                             
tyre_coeffs.pCy1 = P_varGamma(1) ;
tyre_coeffs.pDy1 = P_varGamma(2) ;
tyre_coeffs.pDy2 = P_varGamma(3) ;
tyre_coeffs.pEy1 = P_varGamma(4) ;
tyre_coeffs.pEy2 = P_varGamma(5) ;
tyre_coeffs.pEy3 = P_varGamma(6) ;
tyre_coeffs.pEy4 = P_varGamma(7) ;
tyre_coeffs.pHy1 = P_varGamma(8) ; 
tyre_coeffs.pHy2 = P_varGamma(9) ;
tyre_coeffs.pKy1 = P_varGamma(10) ;
tyre_coeffs.pKy2 = P_varGamma(11) ;
tyre_coeffs.pKy3 = P_varGamma(12) ;
tyre_coeffs.pVy1 = P_varGamma(13) ;
tyre_coeffs.pVy2 = P_varGamma(14) ;

FY0_varGamma_vec = MF96_FY0_vec(zeros_vec, ALPHA_vec, GAMMA_vec, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);

figure('Name','10. Fx0 vs Gamma')
plot(ALPHA_vec,TDataGamma.FY,'o')
hold on
plot(ALPHA_vec,FY0_varGamma_vec,'-')
xlabel('$\alpha$ [-]')
ylabel('$F_{y0}$ [N]')
% Calculate the residuals with the optimal solution found above
res_Fy0_varGamma  = resid_pure_Fy_varGamma(P_varGamma,FY_vec, ALPHA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs);

% R-squared is 
% 1-SSE/SST
% SSE/SST = res_Fx0_nom

% SSE is the sum of squared error,  SST is the sum of squared total
fprintf('R-squared = %6.3f\n',1-res_Fy0_varGamma);


[alpha__y, By, Cy, Dy, Ey, SVy] = MF96_FY0_coeffs(0, 0, GAMMA_vec(3), tyre_coeffs.FZ0, tyre_coeffs);
% 
fprintf('By      = %6.3f\n',By);
fprintf('Cy      = %6.3f\n',Cy);
fprintf('muy      = %6.3f\n',Dy/tyre_coeffs.FZ0);
fprintf('Ey      = %6.3f\n',Ey);
fprintf('SVy     = %6.3f\n',SVy);
fprintf('alpha_y = %6.3f\n',alpha__y);
fprintf('Ky      = %6.3f\n',By*Cy*Dy/tyre_coeffs.FZ0);

% % Longitudinal stiffness
% Kx_vec = zeros(size(load_vec));
% for i = 1:length(load_vec)
%   [kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, 0, load_vec(i), tyre_data);
%   Kx_vec(i) = Bx*Cx*Dx/tyre_data.Fz0;
% end
% 
% figure('Name','Kx vs Fz')
% plot(load_vec,Kx_vec,'o-')

%% Fit coefficients for combined case - $F_z$ = $F_z0$ , $\gamma=0$

% Create tyre dataset for LONG data
vec_samples = 1:1:length(smpl_range_LONG);
tyre_data_LONG = table(); % create empty table

% store raw data in table
tyre_data_LONG.SL =  data_LONG.SL(smpl_range_LONG);
tyre_data_LONG.SA =  data_LONG.SA(smpl_range_LONG)*to_rad;
tyre_data_LONG.FZ =  -data_LONG.FZ(smpl_range_LONG);  % 0.453592  lb/kg
tyre_data_LONG.FX =  data_LONG.FX(smpl_range_LONG);
tyre_data_LONG.FY =  -data_LONG.FY(smpl_range_LONG);
tyre_data_LONG.MZ =  -data_LONG.MZ(smpl_range_LONG);
tyre_data_LONG.IA =  data_LONG.IA(smpl_range_LONG)*to_rad;

% Select intervals
% -12 -6 0 6 12
SL_tol = 2*to_rad;
idx.SL_0    =  0-SL_tol            < tyre_data_LONG.SL & tyre_data_LONG.SL < 0+SL_tol;
idx.SL_12pos =  (12*to_rad-SL_tol)   < tyre_data_LONG.SL & tyre_data_LONG.SL < 12*to_rad+SL_tol;
idx.SL_12neg =  -(12*to_rad+SL_tol)  < tyre_data_LONG.SL & tyre_data_LONG.SL < -12*to_rad+SL_tol;
idx.SL_6pos = (6*to_rad-SL_tol)  < tyre_data_LONG.SL & tyre_data_LONG.SL < 6*to_rad+SL_tol;
idx.SL_6neg = -(6*to_rad+SL_tol) < tyre_data_LONG.SL & tyre_data_LONG.SL < -6*to_rad+SL_tol;
SL_0     = tyre_data_LONG( idx.SL_0, : );
SL_12pos  = tyre_data_LONG( idx.SL_12pos, : );
SL_12neg  = tyre_data_LONG( idx.SL_12neg, : );
SL_6pos  = tyre_data_LONG( idx.SL_6pos, : );
SL_6neg  = tyre_data_LONG( idx.SL_6neg, : );

figure(11)
plot(tyre_data_LONG.SL*to_deg)
hold on
plot(vec_samples(idx.SL_0),   SL_0.SL*to_deg,'.');
plot(vec_samples(idx.SL_12pos),SL_12pos.SL*to_deg,'.');
plot(vec_samples(idx.SL_12neg),SL_12neg.SL*to_deg,'.');
plot(vec_samples(idx.SL_6pos),SL_6pos.SL*to_deg,'.');
plot(vec_samples(idx.SL_6neg),SL_6neg.SL*to_deg,'.');
title('Long slip')
xlabel('Samples [-]')
ylabel('[rad]')

% Intersect LAT data
[TDataComb, ~] = intersect_table_data( tyre_data, GAMMA_0, FZ_700);

% First guess & boundaries
%    [  rBy1     rBy2         rBy3     rCy1        rHy1        rVy1     rVy4      rVy5      rVy6]
P0 = [14.796, 13.1203, -0.00135557, 1.16194, 0.00449364, -0.0115607, 24.6598, -3.65873, -4.61654];
lb = [];
ub = [];

ALPHA_vec = TDataComb.SA;
KAPPA_vec = SL_0.SL;
FY_vec    = TDataComb.FY;
SA_vec = -0.3:0.001:0.3;

[P_comb_fz_nom, ~] = fmincon(@(P)resid_comb_Fy(P,FY_vec,KAPPA_vec,ALPHA_vec,0,FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Update tyre data with new optimal values                             
tyre_coeffs.rBy1 = P_comb_fz_nom(1) ;
tyre_coeffs.rBy2 = P_comb_fz_nom(2) ;
tyre_coeffs.rBy3 = P_comb_fz_nom(3) ;
tyre_coeffs.rCy1 = P_comb_fz_nom(4) ;
tyre_coeffs.rHy1 = P_comb_fz_nom(5) ;
tyre_coeffs.rVy1 = P_comb_fz_nom(6) ;
tyre_coeffs.rVy4 = P_comb_fz_nom(7) ;
tyre_coeffs.rVy5 = P_comb_fz_nom(8) ; 
tyre_coeffs.rVy6 = P_comb_fz_nom(9) ;

FY_comb_fz_nom_vec = MF96_FY_vec(SA_vec , SA_vec , zeros(size(SA_vec)), ...
                              FZ0.*ones(size(SA_vec)),tyre_coeffs);

figure(12)
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
plot(SA_vec,FY_comb_fz_nom_vec,'-','LineWidth',2)
xlabel('$\alpha$ [-]')
ylabel('$F_{Y}$ [N]')

res_Fy=resid_comb_Fy(P_comb_fz_nom,FY_vec,KAPPA_vec,ALPHA_vec,0,tyre_coeffs.FZ0, tyre_coeffs);

fprintf('R-squared_FzNom = %6.3f\n',1-res_Fy);
