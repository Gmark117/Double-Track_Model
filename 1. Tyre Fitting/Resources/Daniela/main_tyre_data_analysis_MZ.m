%% Initialisation
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

data_set = 'Hoosier_B1464run23'; % pure lateral forces
%data_set = 'Hoosier_B1464run30';  % braking/traction (pure log. force) + combined

% tyre geometric data:
% Hoosier	18.0x6.0-10
% 18 diameter in inches
% 6.0 section width in inches
% tread width in inches
diameter = 18*2.56; %
Fz0 = 700;   % [N] nominal load is given
R0  = diameter/2/100; % [m] get from nominal load R0 (m) *** TO BE CHANGED ***

fprintf('Loading dataset ...')
switch data_set
  case 'Hoosier_B1464run23'
  load ([data_set_path, 'Hoosier_B1464run23.mat']); % pure lateral
  cut_start = 27760;
  cut_end   = 54500;
  case 'Hoosier_B1464run30'
  load ([data_set_path, 'Hoosier_B1464run30.mat']); % pure longitudinal
  cut_start = 19028;
  cut_end   = 37643;
  otherwise 
  error('Not found dataset: `%s`\n', data_set) ;
  
end

% select dataset portion
smpl_range = cut_start:cut_end;

fprintf('completed!\n')

%% Plot raw data

figure
tiledlayout(6,1)

ax_list(1) = nexttile; y_range = [min(min(-FZ),0) round(max(-FZ)*1.1)];
plot(-FZ)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Vertical force')
xlabel('Samples [-]')
ylabel('[N]')

ax_list(2) = nexttile; y_range = [min(min(IA),0) round(max(IA)*1.1)];
plot(IA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Camber angle')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(3) = nexttile; y_range = [min(min(SA),0) round(max(SA)*1.1)];
plot(SA)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Side slip')
xlabel('Samples [-]')
ylabel('[deg]')

ax_list(4) = nexttile; y_range = [min(min(SL),0) round(max(SL)*1.1)];
plot(SL)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Longitudinal slip')
xlabel('Samples [-]')
ylabel('[-]')

ax_list(5) = nexttile; y_range = [min(min(P),0) round(max(P)*1.1)];
plot(P)
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre pressure')
xlabel('Samples [-]')
ylabel('[psi]')

ax_list(6) = nexttile;  y_range = [min(min(TSTC),0) round(max(TSTC)*1.1)];
plot(TSTC,'DisplayName','Center')
hold on
plot(TSTI,'DisplayName','Internal')
plot(TSTO,'DisplayName','Outboard')
hold on
plot([cut_start cut_start],y_range,'--r')
plot([cut_end cut_end],y_range,'--r')
title('Tyre temperatures')
xlabel('Samples [-]')
ylabel('[degC]')

linkaxes(ax_list,'x')


%plot(SA,FY)

%% Select some specific data
% Cut crappy data and select only 12 psi data

vec_samples = 1:1:length(smpl_range);
tyre_data = table(); % create empty table
% store raw data in table
tyre_data.SL =  SL(smpl_range);
tyre_data.SA =  SA(smpl_range)*to_rad;
tyre_data.FZ =  -FZ(smpl_range);  % 0.453592  lb/kg
tyre_data.FX =  FX(smpl_range);
tyre_data.FY =  -FY(smpl_range);
tyre_data.MZ =  -MZ(smpl_range);
tyre_data.IA =  IA(smpl_range)*to_rad;

% Extract points at constant inclination angle
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
% Test data done at: 
%  - 50lbf  ( 50*0.453592*9.81 =  223N )
%  - 150lbf (150*0.453592*9.81 =  667N )
%  - 200lbf (200*0.453592*9.81 =  890N )
%  - 250lbf (250*0.453592*9.81 = 1120N )

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
SA_tol = 0.5*to_rad;
idx.SA_0    =  0-SA_tol          < tyre_data.SA & tyre_data.SA < 0+SA_tol;
idx.SA_3neg = -(3*to_rad+SA_tol) < tyre_data.SA & tyre_data.SA < -3*to_rad+SA_tol;
idx.SA_6neg = -(6*to_rad+SA_tol) < tyre_data.SA & tyre_data.SA < -6*to_rad+SA_tol;
SA_0     = tyre_data( idx.SA_0, : );
SA_3neg  = tyre_data( idx.SA_3neg, : );
SA_6neg  = tyre_data( idx.SA_6neg, : );

figure()
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
hold on
plot(vec_samples(idx.SA_0),   SA_0.SA*to_deg,'.');
plot(vec_samples(idx.SA_3neg),SA_3neg.SA*to_deg,'.');
plot(vec_samples(idx.SA_6neg),SA_6neg.SA*to_deg,'.');
title('Slide slip')
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

figure('Name','Selected-data')
plot_selected_data(TData0);

%% FITTING 
% initialise tyre data
tyre_coeffs = initialise_tyre_data(R0, Fz0);

%% Fitting with Fz=Fz_nom= 700N and camber=0  alpha = 0 VX= 10
% ------------------
% self aligning torque

% Fit the coeffs {qHz1, qBz1, qCz1, qDz1, qEz1, qEz4, qBz9, qBz10, qDz6}
FZ0 = mean(TData0.FZ);

zeros_vec = zeros(size(TData0.SA));
ones_vec  = ones(size(TData0.SA));

MZ0_guess = MF96_MZ0_vec(zeros_vec, TData0.SA, zeros_vec, tyre_coeffs.FZ0*ones_vec, tyre_coeffs);

% check guess 
figure()
plot(TData0.SA,TData0.MZ,'.')
hold on
plot(TData0.SA,MZ0_guess,'-')

% Plot raw data and initial guess
% figure()
% plot(TDataSub.KAPPA,TDataSub.FX,'o')
% hold on
% plot(TDataSub.KAPPA,FX0_guess,'x')

% Guess values for parameters to be optimised
%    [qBz1 qBz9 qBz10 qCz1 qDz1 qDz6 qEz1 qEz4 qHz1] 
% P0 = [ 10, -500, -500,  2, 0.1, 0.1, 1.5,   0,    0]; 
P0 = [ 10, 50, 0,  1, 0.1, 0, -1.5,   0,    0];
% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
% [qBz1  qBz9  qBz10  qCz1  qDz1  qDz6  qEz1  qEz4  qHz1]
% lb = [  9,  -501, -100,   1,    0,    0,    1, -0.1,    -1];
% ub = [ 20,  -499,  100,   3,  0.2,    1,    2,    1,     1];
% lb = [-10, -500,  -1000, -10, -10, -10, -10, -10, -10];
% ub = [ 50,  500,   1000,  10,  10,  10,  10,  10,  10];
% lb = [  90,  -501, -1000,   -100, -10, -100, -100, 0, 0];
% ub = [ 200,  -499,  1000,   300,  20,    100,    100, 0, 0];
lb = [  10,  20, -1,   10, -10,   -10,   -20, -10,    -10];
ub = [ 13,  40,  1,   20,  1,    10,    10,  10,    10];


ALPHA_vec = TData0.SA;
MZ_vec    = TData0.MZ;

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

[P_fz_nom, ~] = fmincon(@(P)resid_pure_Mz(P,MZ_vec, ALPHA_vec,0,tyre_coeffs.FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

%    [qBz1  qBz9  qBz10  qCz1  qDz1  qDz6  qEz1  qEz4  qHz1]
% Update tyre data with new optimal values  
tyre_coeffs.qBz1  = P_fz_nom(1) ;
tyre_coeffs.qBz9  = P_fz_nom(2) ; 
tyre_coeffs.qBz10 = P_fz_nom(3) ;
tyre_coeffs.qCz1  = P_fz_nom(4) ;
tyre_coeffs.qDz1  = P_fz_nom(5) ;
tyre_coeffs.qDz6  = P_fz_nom(6) ;
tyre_coeffs.qEz1  = P_fz_nom(7) ;
tyre_coeffs.qEz4  = P_fz_nom(8) ;
tyre_coeffs.qHz1  = P_fz_nom(9) ;

MZ0_fz_nom_vec = MF96_MZ0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), ...
                              FZ0.*ones(size(SA_vec)),tyre_coeffs);

figure('Name','Mz0(Fz0)')
plot(TData0.SA,TData0.MZ,'o')
hold on
%plot(TDataSub.KAPPA,FX0_fz_nom_vec,'-')
plot(SA_vec,MZ0_fz_nom_vec,'-','LineWidth',2)
xlabel('$\alpha$ [-]')
ylabel('$M_{z0}$ [N]')

res_MZ0= resid_pure_Mz(P,MZ_vec, ALPHA_vec,0,tyre_coeffs.FZ0, tyre_coeffs);
fprintf('R-squared = %6.3f\n',1-res_MZ0);

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
MZ0_guess = MF96_MZ0_vec(zeros_vec, TDataDFz.SA, zeros_vec, FZ0*ones_vec, tyre_coeffs);

% check guess 
figure(6)
plot(TDataDFz.SA,TDataDFz.MZ,'.')
hold on
plot(TDataDFz.SA,MZ0_guess,'-')
% plot(TDataDFz.SA,FY0_fz_nom_vec,'-','LineWidth',2)

% Plot raw data and initial guess
% figure()
% plot(TDataSub.ALPHA,TDataSub.FY,'o')
% hold on
% plot(TDataSub.ALPHA,FY0_guess,'y')

% Guess values for parameters to be optimised
%    [qBz2  qBz3   qDz2   qDz7    qEz2    qEz3    qHz2]
P0 = [  -2, -0.5, 0.005, 0.01,  -0.005, -0.001, 0.0001];


% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%    []
lb = [];
ub = [];


ALPHA_vec = TDataDFz.SA;
MZ_vec    = TDataDFz.MZ;
FZ_vec    = TDataDFz.FZ;

% check guess
SA_vec = -0.3:0.001:0.3;
MZ0_dfz_vec = MF96_MZ0_vec(zeros(size(SA_vec)), SA_vec, zeros(size(SA_vec)), ...
                           TDataDFz.FZ,tyre_coeffs);
% 
% figure
% plot(KAPPA_vec,FX_vec,'.')
% hold on
% plot(SL_vec,FX0_dfz_vec,'.')


% LSM_pure_Fx returns the residual, so minimize the residual varying X. It
% is an unconstrained minimization problem 

[P_dfz,~,exitflag] = fmincon(@(P)resid_pure_Mz_varFz(P,MZ_vec, ALPHA_vec,0,FZ_vec, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

disp(exitflag)
% Change tyre data with new optimal values 
tyre_coeffs.qBz2 = P_dfz(1) ; % 1
tyre_coeffs.qBz3 = P_dfz(2) ;
tyre_coeffs.qDz2 = P_dfz(3) ;
tyre_coeffs.qDz7 = P_dfz(4) ;
tyre_coeffs.qEz2 = P_dfz(5) ;
tyre_coeffs.qEz3 = P_dfz(6) ;
tyre_coeffs.qHz2 = P_dfz(7) ; 


res_MZ0_dfz_vec = resid_pure_Mz_varFz(P_dfz,MZ_vec,SA_vec,0 , FZ_vec,tyre_coeffs);

tmp_zeros = zeros(size(SA_vec));
tmp_ones = ones(size(SA_vec));


MZ0_fz_var_vec1 = MF96_MZ0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_220.FZ)*tmp_ones,tyre_coeffs);
MZ0_fz_var_vec2 = MF96_MZ0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_440.FZ)*tmp_ones,tyre_coeffs);
MZ0_fz_var_vec3 = MF96_MZ0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_700.FZ)*tmp_ones,tyre_coeffs);
MZ0_fz_var_vec4 = MF96_MZ0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_900.FZ)*tmp_ones,tyre_coeffs);
MZ0_fz_var_vec5 = MF96_MZ0_vec(tmp_zeros, SA_vec ,tmp_zeros, mean(FZ_1120.FZ)*tmp_ones,tyre_coeffs);


figure('Name','7. Mz0(Fz0)')
plot(TDataDFz.SA,TDataDFz.MZ,'o')
hold on
%plot(TDataSub.ALPHA,FY0_fz_nom_vec,'-')
%plot(SA_vec,FY0_dfz_vec,'-','LineWidth',2)
plot(SA_vec,MZ0_fz_var_vec1,'-','LineWidth',2)
plot(SA_vec,MZ0_fz_var_vec2,'-','LineWidth',2)
plot(SA_vec,MZ0_fz_var_vec3,'-','LineWidth',2)
plot(SA_vec,MZ0_fz_var_vec4,'-','LineWidth',2)
plot(SA_vec,MZ0_fz_var_vec5,'-','LineWidth',2)

xlabel('$\alpha$ [-]')
ylabel('$F_{y0}$ [N]')

%% Fit coefficient with variable camber

% extract data with variable load
[TDataGamma, ~] = intersect_table_data( tyre_data, FZ_700 );

% Fit the coeffs { pDx3}

% Guess values for parameters to be optimised
%    [ b4   b4   d3   d4   d8   d9   e5   h3   h4]
P0 = [  -0.01,-0.01,0.01,0.01,0.5,10,0.01,5,-1]; 

% NOTE: many local minima => limits on parameters are fundamentals
% Limits for parameters to be optimised
% 1< pCx1 < 2 
% 0< pEx1 < 1 
%    [ b4   b4   d3   d4   d8   d9   e5   h3   h4]
lb = [];
ub = [];


zeros_vec = zeros(size(TDataGamma.SA));
ones_vec  = ones(size(TDataGamma.SA));

ALPHA_vec = TDataGamma.SA;
GAMMA_vec = TDataGamma.IA; 
MZ_vec    = TDataGamma.MZ;
FZ_vec    = TDataGamma.FZ;

figure(9)
plot(ALPHA_vec,MZ_vec);


% LSM_pure_Fy returns the residual, so minimize the residual varying Y. It
% is an unconstrained minimization problem 
[P_varGamma,fval,exitflag] = fmincon(@(P)resid_pure_Mz_varGamma(P,MZ_vec, ALPHA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs),...
                               P0,[],[],[],[],lb,ub);

% Change tyre data with new optimal values     
%    [ b4   b5   d3   d4   d8   d9   e5   h3   h4]                        
tyre_coeffs.qBz4 = P_varGamma(1) ; % 1
tyre_coeffs.qBz5 = P_varGamma(2) ;
tyre_coeffs.qDz3 = P_varGamma(3) ;
tyre_coeffs.qDz4 = P_varGamma(4) ;
tyre_coeffs.qDz8 = P_varGamma(5) ;
tyre_coeffs.qDz9 = P_varGamma(6) ;
tyre_coeffs.qEz5 = P_varGamma(7) ;
tyre_coeffs.qHz3 = P_varGamma(8) ; 
tyre_coeffs.qHz4 = P_varGamma(9) ; 

MZ0_varGamma_vec = MF96_MZ0_vec(zeros_vec, ALPHA_vec, GAMMA_vec, tyre_coeffs.FZ0*ones_vec,tyre_coeffs);

figure('Name','10. Mz0 vs Gamma')
plot(ALPHA_vec,TDataGamma.MZ,'o')
hold on
plot(ALPHA_vec,MZ0_varGamma_vec,'-')
xlabel('$\alpha$ [-]')
ylabel('$M_{z0}$ [N]')
% Calculate the residuals with the optimal solution found above
res_Mz0_varGamma  = resid_pure_Mz_varGamma(P_varGamma,MZ_vec, ALPHA_vec,GAMMA_vec,tyre_coeffs.FZ0, tyre_coeffs);

% R-squared is 
% 1-SSE/SST
% SSE/SST = res_Fx0_nom

% SSE is the sum of squared error,  SST is the sum of squared total
fprintf('R-squared = %6.3f\n',1-res_Mz0_varGamma);


[alpha__t, alpha__r, Bt, Ct, Dt, Et, Br, Dr] = MF96_MZ0_coeffs(0, 0, GAMMA_vec(3), tyre_coeffs.FZ0, tyre_coeffs);
% 
fprintf('Bt      = %6.3f\n',Bt);
fprintf('Ct      = %6.3f\n',Ct);
fprintf('Dt      = %6.3f\n',Dt);
fprintf('Et      = %6.3f\n',Et);
fprintf('alpha_t = %6.3f\n',alpha__t);
fprintf('Br      = %6.3f\n',Br);
fprintf('Dr      = %6.3f\n',Dr);
fprintf('alpha_r = %6.3f\n',alpha__r);

% % Longitudinal stiffness
% Kx_vec = zeros(size(load_vec));
% for i = 1:length(load_vec)
%   [kappa__x, Bx, Cx, Dx, Ex, SVx] = MF96_FX0_coeffs(0, 0, 0, load_vec(i), tyre_data);
%   Kx_vec(i) = Bx*Cx*Dx/tyre_data.Fz0;
% end
% 
% figure('Name','Kx vs Fz')
% plot(load_vec,Kx_vec,'o-')
