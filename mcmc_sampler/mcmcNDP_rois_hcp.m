%for nDP process
function [zeta, lP_record, xi, saveid] = mcmcNDP_rois_hcp(ith_run,idx_roi,In_COMP2,COMP1,ALPHA_BETA_PRIOR)
% function for perform MCMC NDP sampling using the HCP data

% input: ith_run - index for this run ?for saving purpose?
%        idx_roi - index for rois; in the range of 1-45
%        In_COMP2 - indicator for including the translation component
%        COMP1 - name of component 1, can have two values,'fpca' or 'trans'
%        ALPHA_BETA_PRIOR - indicator to set priors for alpha and beta
    
% output: zeta: MCMC samples of membership indicator for the each individual
%         lP_record: log probability record
%         xi: MCMC samples of membership indicator for each fiber curve in each subject


%parameters
Npoint = 100;
Kcluster = 1;

%parameter for subsampling
SubSamp =1;
Nsubsample = 600;

%load data
load(sprintf('./data/selected_hcp_rois_fpcaresults.mat'))

ssubid = [6,3,4,5,14,17,18,16]; %select a subset of subjects

for b=1:8
    n = ssubid(b);
    
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
    inputdata{b}.Orgfiber = newfibers;
    inputdata{b}.Rots = Rots(:,:,AllInD);
    inputdata{b}.Trans = Trans(:,AllInD);
    
    N = length(AllInD);
    fpca_coeff = [];
    for i=1:N
        id = AllInD(i);
        tmp_coeff = PCA_coeffs{id};
        coeff_vect = [tmp_coeff{1}(1),tmp_coeff{2}(1),tmp_coeff{3}(1)]';
        fpca_coeff(:,i) = coeff_vect;
        new_PCA_coeffs{i} = tmp_coeff;
    end
    inputdata{b}.fpca_coeff = fpca_coeff;
    
    for i=1:N
        R = Rots(:,:,i);
        trR = trace(R);
        theta(i) = acos((trR-1)/2);
        tempM = (R - R')*theta(i)/(2*sin(theta(i)));
        mapedR(:,i) = [-R(1,2),R(1,3),R(3,2)]';
    end
    inputdata{b}.mapedR = mapedR;
    inputdata{b}.Rots = mapedR;
    
    %reconstruction
    %final_reconstl = roifibers_reconstruction_adaptive(labels(AllInD),Trans(:,AllInD),Rots(:,:,AllInD),new_PCA_coeffs,means_sl,pca_basis);
    %inputdata{n}.Reconfiber = newfibers;
    clear new_PCA_coeffs;
    
    % remove mean from the trans
    inputdata{b}.Trans = inputdata{b}.Trans - repmat(mean(inputdata{b}.Trans,2),1,size(inputdata{b}.Trans,2));
end

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



option.INDICATOR_COMP3=0;
option.INDICATOR_COMP2=In_COMP2;
option.COMP1 = COMP1;
option.NIW = 1;
option.GRAPHICS=0;
option.ALPHA_BETA_PRIOR = ALPHA_BETA_PRIOR;
option.FIRST_LS = 1; % first label-swithing
option.SECOND_LS = 1; %second label-swithing;
option.MCMCdiag = 0;

[zeta,lP_record,xi] = samplerNDP(inputdata,8,10,8,10,5000,option);

eval(sprintf('save mcmcresult_hcp_roi%d_%d zeta lP_record xi option',idx_roi,ith_run));

saveid = 1;
