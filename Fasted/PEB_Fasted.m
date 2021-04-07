
%% SECOND LEVEL SPDCM - FOOD CHOICE CIRCUIT
GCM = {'DCM_swr_31003f.mat';'DCM_swr_31004f.mat';'DCM_swr_31011f.mat';'DCM_swr_31015f.mat';'DCM_swr_31016f.mat';'DCM_swr_31017f.mat';'DCM_swr_31019f.mat';'DCM_swr_31020f.mat';...
       'DCM_swr_31022f.mat';'DCM_swr_31023f.mat';'DCM_swr_31026f.mat';'DCM_swr_31027f.mat';'DCM_swr_31031f.mat';'DCM_swr_31032f.mat';'DCM_swr_31033f.mat';'DCM_swr_31036f.mat';...
       'DCM_swr_31037f.mat';'DCM_swr_31042f.mat';'DCM_swr_31047f.mat';'DCM_swr_32006f.mat';'DCM_swr_32007f.mat';'DCM_swr_32009f.mat';'DCM_swr_32012f.mat';'DCM_swr_32013f.mat';...
       'DCM_swr_32014f.mat';'DCM_swr_32018f.mat';'DCM_swr_32021f.mat';'DCM_swr_32024f.mat';'DCM_swr_32025f.mat';'DCM_swr_32028f.mat';'DCM_swr_32029f.mat';'DCM_swr_32030f.mat';...
       'DCM_swr_32034f.mat';'DCM_swr_32038f.mat';'DCM_swr_32041f.mat';'DCM_swr_32043f.mat';'DCM_swr_32044f.mat';'DCM_swr_32046f.mat'}

addpath(genpath('/Users/katharinavoigt/Documents/MATLAB/spm12'))
path = pwd;


M = struct();
M.alpha = 1;
M.beta  = 16;
M.hE    = 0;
M.hC    = 1/16;
M.Q     = 'single';

N = length(GCM); % no. of subjects

% Specify design matrix for N subjects. It should start with a constant column
cd('/Users/katharinavoigt/Library/Mobile Documents/com~apple~CloudDocs/Research/POSTDOC PAPERS/05.Hypothalamus/NEWANALYSIS/Results_N38/BehaviouralData')
load('./AllBehaviourals.mat');
cd(path);

% Define design matric
M.X = [ones(N,1)];

% Choose field
field = {'A'};

% Estimate model
PEB     = spm_dcm_peb(GCM,M,field);

BMA = spm_dcm_peb_bmc(PEB);

spm_dcm_peb_review(BMA,GCM);

save('Model1_Intercept');

% % Perform leave-one-out cross validation (GCM,M,field are as before)
spm_dcm_loo(GCM,M,field);

save('Model1_Intercept');

%% Plotting stuff after
spm_dcm_fmri_check();
spm_dcm_fmri_csd_results();

% first level plots / checks
% 2nd level plots
n = 4; % number of regions
Ep = PEB.Ep(:,2); % using the second covariate
Vp = diag(PEB.Cp(1:n*n,1:n*n));

figure;
spm_plot_ci(Ep,Vp);

% If you want to calculate the postrerior probabilities then you can simply type

Pp = 1 - spm_Ncdf(0,abs(Ep),Vp); % Ep and Vp from above

%which will give you probabilities for exampe 0.96 is equal to 96% probability (or confidence) that that connection was larger than zero.

%If you want to calculate the 90% credible intervals then type

ci = spm_invNcdf(0.95,Ep,Vp);
upper = Ep + ci;
lower = Ep - ci;


