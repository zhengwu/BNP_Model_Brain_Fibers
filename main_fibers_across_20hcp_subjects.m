% for modeling fiber curves between two regions of interet in multiple subjects
% this main program will save the posterior results into a .mat file
% here we use 20 subjects' data from the HCP dataset

clear all;
close all;

% add necessary paths
addpath('./data');
addpath('./utilities');
addpath('./curve_processing_functions')
addpath('./mcmc_sampler')
addpath(genpath('./MCMCDiag'))

if(exist('mydistributions','dir')==7)
    %do nothing
else
    %unzip the file
     unzip('mydistributions.zip');
     addpath(genpath('./mydistributions'))
end

%%%%%%%%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%
load(sprintf('selected_rois.mat'));


%load data (big data, can't host over the gitbuh repository, need to download first)
if(exist('./data/selected_hcp_rois_fpcaresults.mat','file')==2)
    %do nothing
else
    %download the data
    outfilename = websave('./data/selected_hcp_rois_fpcaresults.mat','https://www.dropbox.com/s/qdc4qv21tv9oboj/selected_hcp_rois_fpcaresults.mat?dl=1');
end



%%%%%%%%%%%%%%%%%%%%%% set parameters %%%%%%%%%%%%%%%
ith_run = 1; % index of the run (for saving the final results)
idx_roi = 45; % index of pair of regions (ranging from 1 to 45)
roia = array_a(idx_roi);
roib = array_b(idx_roi);
display(sprintf('running the posterior sample for roia_%d roib_%d',roia,roib));



%%%%%%%%%%%%%%%%%%%%%% perform posterior sampling %%%%%%%%%%%%%%%%%%%%%
Indicator_COMP2 = 0; % indicator for combining the second component
% 1 - include the second component (the default second component is translation component) 
% 0 - only use COMP1 for doing the posterior inference

COMP1 = 'fpca'; %  'fpca', choose function pca coefficents (the shape component) as the first component 
                % 'trans', choose the translation component as the first component                
Indicator_ALPHA_BETA_PRIOR = 1; % set priors for parameters alpha and beta in NDP;

% seed for reproducibility
% load('random_seed.mat');
% rng(d);
[zeta, lP_record, xi, saveid] = mcmcNDP_rois_hcp(ith_run,idx_roi,0,COMP1,Indicator_ALPHA_BETA_PRIOR); 

%%%%%%%%%%%%%%%%%%%%%% process and analysis the final result %%%%%%%%%%%%%%%%%%%%%
Nsample = size(zeta,1);
N = 20;
burnin = 500;
finalB = zeros(N,N);
for i=burnin:Nsample
    tempB = calculate_B(squeeze(zeta(i,:)),N);
    finalB = finalB + tempB;
    psoteriorK(i-burnin+1) = length(unique(zeta(i,:)));
end

%plot figure
figure(1);clf; imagesc(finalB);
set(gca,'fontsize',20);
figure(2),clf; hist(psoteriorK)

mode(psoteriorK)
[c,cn,thrd] = BtoCluster_zw(finalB,N,mode(psoteriorK));
[AR,RI,MI,HI]=RandIndex(c,[ones(1,10),2*ones(1,10)]);