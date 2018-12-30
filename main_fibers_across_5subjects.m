% for modeling fiber curves between two regions of interet in multiple subjects
% this main program will save the posterior results into a .mat file
% here we use 5 subjects' data from the Sherbrooke test retest dataset

clear all;
close all;

%%%%%%%%%%%%%%%% load and preprocess data %%%%%%%
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

%parameters
Npoint = 100;
Kcluster = 1;

%parameter for subsampling
SubSamp = 1;
Nsubsample = 400;

% select which connection to study (refer to the paper for more details on
% each connection)
load(sprintf('selected_rois.mat'));
idx_roi = 2; % index of pair of regions (1-45)
roia = array_a(idx_roi);
roib = array_b(idx_roi);

%load data
load(sprintf('./data/selected_rois_fpcaresults_nofibers.mat'))

AllRotation = [];
for n=1:15
    Fpcaresults = selected_rois_fpcaresults{idx_roi,n};
    labels = Fpcaresults.labels;
    Trans = Fpcaresults.Trans;
    Rots = Fpcaresults.Rots;
    PCA_coeffs = Fpcaresults.PCA_coeffs;
    
    Nfiber = length(labels);
    AllInD = 1:Nfiber;
    if(SubSamp == 1)
        sampledidx = datasample(AllInD, min(Nsubsample,Nfiber));
        AllInD = sampledidx;
    end
    
    %newfibers = normal_fibers(:,:,AllInD);
    newfibers=[];
    inputdata{n}.Orgfiber = newfibers;
    inputdata{n}.Rots = Rots(:,:,AllInD);
    inputdata{n}.Trans = Trans(:,AllInD);
    
    AllRotation = cat(3, AllRotation, inputdata{n}.Rots);
    
    N = length(AllInD);
    fpca_coeff = [];
    for i=1:N
        id = AllInD(i);
        tmp_coeff = PCA_coeffs{id};
        coeff_vect = [tmp_coeff{1}(1),tmp_coeff{2}(1),tmp_coeff{3}(1)]';
        fpca_coeff(:,i) = coeff_vect;
        new_PCA_coeffs{i} = tmp_coeff;
    end;
    inputdata{n}.fpca_coeff = fpca_coeff;
    
    %         for i=1:N
    %             R = Rots(:,:,i);
    %             trR = trace(R);
    %             theta(i) = acos((trR-1)/2);
    %             tempM = (R - R')*theta(i)/(2*sin(theta(i)));
    %             mapedR(:,i) = [-R(1,2),R(1,3),R(3,2)]';
    %         end
    %         inputdata{n}.mapedR = mapedR;
    
    clear new_PCA_coeffs;
    
    % remove mean from the trans
    inputdata{n}.Trans = inputdata{n}.Trans - repmat(mean(inputdata{n}.Trans,2),1,size(inputdata{n}.Trans,2));
end

% preprocess rotation;
Omean = Mean_Rotation(AllRotation);

for n=1:15
    Rots = inputdata{n}.Rots;
    for i=1:size(Rots,3)
        newRots(:,:,i) = Omean'*Rots(:,:,i);
        
        R = newRots(:,:,i);
        trR = trace(R);
        theta(i) = acos((trR-1)/2);
        tempM = (R - R')*theta(i)/(2*sin(theta(i)));
        mapedR(:,i) = [-R(1,2),R(1,3),R(3,2)]';
    end
    inputdata{n}.Rots = mapedR;
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
option.INDICATOR_COMP2=0;
option.COMP1 = 'fpca';
option.NIW = 1;
option.GRAPHICS=0;
option.ALPHA_BETA_PRIOR = 1;
option.FIRST_LS = 1; %first label-swithing;
option.SECOND_LS = 1; %second label-swithing;
option.MCMCdiag = 0;

[zeta,lP_record,xi] = samplerNDP(inputdata,9,15,9,15,5000,option);

%idx_roi = 2;
%ith_run = 1;
%eval(sprintf('save mcmcresult_roi_5sub_rotation_%d_%d zeta lP_record xi option',idx_roi,ith_run));

%%%%%%%%%%%%%%%%%%%%%% process and analysis the final result %%%%%%%%%%%%%%%%%%%%%
Nsample = size(zeta,1);
N = 15;
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
[AR,RI,MI,HI]=RandIndex(c,[1,1,1,2,2,2,3,3,3,4,4,4,5,5,5]);