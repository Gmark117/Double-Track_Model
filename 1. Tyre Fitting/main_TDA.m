%% Initialisation
clc
clearvars
close all

addpath('dataset/')
addpath('tyre_lib/')
addpath('plots/')

% Set LaTeX as default interpreter for axis labels, ticks and legends
set(0,'defaulttextinterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultAxesFontSize',  16)
set(0,'DefaultLegendFontSize',16)

to_rad = pi/180;
to_deg = 180/pi;

%% Data Loading
dataset_LAT  = 'Hoosier_B1464run23';  % Pure Lateral
dataset_LONG = 'Hoosier_B1464run30';  % Pure Longitudinal
data_LAT     = load(dataset_LAT);
data_LONG    = load(dataset_LONG);

% Cut data
cut_start_LAT     = 27760;
cut_end_LAT       = 54500;
cut_start_LONG    = 19028;
cut_end_LONG      = 37643;
sample_range_LAT  = cut_start_LAT:cut_end_LAT;
sample_range_LONG = cut_start_LONG:cut_end_LONG;

% Store data
vec_samples_LAT  = 1:1:length(sample_range_LAT);
vec_samples_LONG = 1:1:length(sample_range_LONG);
tyre_data_LAT    = table();
tyre_data_LONG   = table();

tyre_data_LAT.SL  =  data_LAT.SL(sample_range_LAT);
tyre_data_LAT.SA  =  data_LAT.SA(sample_range_LAT)*to_rad;
tyre_data_LAT.FZ  =  -data_LAT.FZ(sample_range_LAT);
tyre_data_LAT.FX  =  data_LAT.FX(sample_range_LAT);
tyre_data_LAT.FY  =  -data_LAT.FY(sample_range_LAT);
tyre_data_LAT.MZ  =  -data_LAT.MZ(sample_range_LAT);
tyre_data_LAT.IA  =  data_LAT.IA(sample_range_LAT)*to_rad;

tyre_data_LONG.SL =  data_LONG.SL(sample_range_LONG);
tyre_data_LONG.SA =  data_LONG.SA(sample_range_LONG)*to_rad;
tyre_data_LONG.FZ =  -data_LONG.FZ(sample_range_LONG);
tyre_data_LONG.FX =  data_LONG.FX(sample_range_LONG);
tyre_data_LONG.FY =  data_LONG.FY(sample_range_LONG);
tyre_data_LONG.MZ =  data_LONG.MZ(sample_range_LONG);
tyre_data_LONG.IA =  data_LONG.IA(sample_range_LONG)*to_rad;

%% Raw Data Plotting
plot_raw_data(data_LAT, cut_start_LAT, cut_end_LAT, 'Lateral Raw Data')
plot_raw_data(data_LONG, cut_start_LONG, cut_end_LONG, 'Longitudinal Raw Data')

%% Subtables Creation
Tsubs_LAT = generate_subtables(tyre_data_LAT, vec_samples_LAT, 'Lateral Subtables');
Tsubs_LONG = generate_subtables(tyre_data_LONG, vec_samples_LONG, 'Longitudinal Subtables');

%% Data Initialisation
tyre_coeffs = initialise_tyre_data();

diameter = 18*2.56;
tyre_coeffs.R0 = diameter/2/100;

%% Fitting
fitting_LONG;
fitting_LAT;
fitting_ALIGN;
