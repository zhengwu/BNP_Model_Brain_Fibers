%for modeling fiber curves between two regions of interet in multiple subjects
% here we use 3 subjects' data from the Sherbrooke test retest dataset


clear all;
close all;

%%%%%%%%%%%%%%%% load data %%%%%%%
Nsamples = 5000;
roia = 16;
roib = 55;

addpath(genpath('./utilities'));
addpath(genpath('./curve_processing_functions'));
addpath(genpath('./data'))
addpath(genpath('./MCMCDiag'))
addpath('./mcmc_sampler')

if(exist('mydistributions','dir')==7)
    %do nothing
else
    %unzip the file
     unzip('mydistributions.zip');
     addpath(genpath('./mydistributions'))
end


%load data
load(sprintf('./data/roi16_55_rotation.mat'))

% %%%%%%%%%%% get the data ready %%%%%%%%%%
idx = 0;
AllRotation = [];
for idsub = 1:3
    for iscan = 1:3
        
        idx = idx + 1;
        newfibers = inputdata{idx}.Orgfiber;
        Nfiber = size(newfibers,3); % # of fibers
        Rots = inputdata{idx}.Rots;
        Trans = inputdata{idx}.Trans;
        fpca_coeff = inputdata{idx}.fpca_coeff;
        final_reconstl = inputdata{idx}.Reconfiber;
        
        AllRotation = cat(3, AllRotation, Rots);
        
    end
end

%%%%%%%%%%%%%%%%%%% process rotation %%%%%%%%%%%%%%%%%%%%
% preprocess rotation;
Omean = Mean_Rotation(AllRotation);

for idx=1:9
    Rots = inputdata{idx}.Rots;
    for i=1:size(Rots,3)
        newRots(:,:,i) = Omean'*Rots(:,:,i);
        
        R = newRots(:,:,i);
        trR = trace(R);
        theta(i) = acos((trR-1)/2);
        tempM = (R - R')*theta(i)/(2*sin(theta(i)));
        mapedR(:,i) = [-R(1,2),R(1,3),R(3,2)]';
    end
    inputdata{idx}.Rots = mapedR;
    clear mapedR;
end

% %%%%%%%%%%% perform posterior sampling using NDP %%%%%%%%%%
% option here controls most of the parameters inside NDP
% INDICATOR_COMP2 - whether to include component 2
% INDICATOR_COMP3 - wheter to include component 3
% COMP1 - string, should be 'fpca' or 'trans' in with the inputdata
%         'fpca' indicates that the first component is the functional pca coefficients
%         'trans' indicates that the first component is the translation
% NIW - indicator for using normal-inverse-Wishart distribution as the
% prior for theta
% GRAPHICS - indicator for displaying the MCMC chain while sampling
% ALPHA_BETA_PRIOR - indicator for setting priors for alpha and beta in NDP.
% FIRST_LS & SECOND_LS - indicator for using the label-switching technique (by Papaspiliopoulos & Roberts, 2008) to improve the MCMC sampling efficiency 
% MCMCdiag - indicator for MCMC diagnosis

option.INDICATOR_COMP3=0;
option.INDICATOR_COMP2=0; % we only use the first component 
option.COMP1 = 'fpca'; % set the first component to be the functional pca coefficents


option.NIW = 1;
option.GRAPHICS=0;
option.ALPHA_BETA_PRIOR = 1;
option.FIRST_LS = 1; % first label-swithing
option.SECOND_LS = 0; %second label-swithing;
option.MCMCdiag = 0;

% seed for reproducibility
% load('random_seed.mat');
% rng(d);

[zeta,lP_record,xi] = samplerNDP(inputdata,9,15,9,15,Nsamples,option);

%idx_roi = 1655;
%ith_run = 1;
%eval(sprintf('save mcmcresult_roi_3sub_rotation_%d_%d zeta lP_record xi option',idx_roi,ith_run));