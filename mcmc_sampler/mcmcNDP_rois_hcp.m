%for nDP process
function [zeta, lP_record, xi, saveid] = mcmcNDP_rois(ith_run,idx_roi,In_COMP2,COMP1,ALPHA_BETA_PRIOR)
    addpath(genpath('./distributions'));
    addpath(genpath('./utilities'));
    addpath(genpath('./curve_processing_functions'));
    
    %parameters
    Npoint = 100;
    Kcluster = 1;
    
    %perfor the subsampling
    SubSamp =1;
    Nsubsample = 200;

   
    %load data
    load(sprintf('./data/selected_rois_fpcaresults_nofibers.mat'))
    %load(sprintf('../data/selected_rois_fpcaresults_rotation.mat'))
    
    for n=1:15
        
        Fpcaresults = selected_rois_fpcaresults{idx_roi,n};
    
        labels = Fpcaresults.labels;
        Trans = Fpcaresults.Trans;
        Rots = Fpcaresults.Rots;
        PCA_coeffs = Fpcaresults.PCA_coeffs;
        %normal_fibers = Fpcaresults.normal_fibers;
        
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
        
        for i=1:N
            R = Rots(:,:,i);
            trR = trace(R);
            theta(i) = acos((trR-1)/2);
            tempM = (R - R')*theta(i)/(2*sin(theta(i)));
            mapedR(:,i) = [-R(1,2),R(1,3),R(3,2)]';
        end
        inputdata{n}.Rots = mapedR;
        
          %reconstruction
          %final_reconstl = roifibers_reconstruction_adaptive(labels(AllInD),Trans(:,AllInD),Rots(:,:,AllInD),new_PCA_coeffs,means_sl,pca_basis);
          %inputdata{n}.Reconfiber = newfibers;
          clear new_PCA_coeffs;
    
          % remove mean from the trans
          inputdata{n}.Trans = inputdata{n}.Trans - repmat(mean(inputdata{n}.Trans,2),1,size(inputdata{n}.Trans,2));
    end
   
    
    option.INDICATOR_COMP3=0;
    option.INDICATOR_COMP2=In_COMP2;
    option.COMP1 = COMP1;
    option.NIW = 1;
    option.GRAPHICS=0;
    option.ALPHA_BETA_PRIOR = ALPHA_BETA_PRIOR;
    option.FIRST_LS = 1; %first label-swithing;
    option.SECOND_LS = 1; %second label-swithing;

    [zeta,lP_record,xi] = samplerNDP(inputdata,15,20,15,20,5000,option);
    %[zeta,lP_record,xi] = samplerNDP_memlimited(inputdata,15,20,15,20,5000,option);
    %[zeta,lP_record,xi] = samplerNDP_distribution(inputdata,9,15,9,15,5000,option);
    
    eval(sprintf('save noRot_mcmcresult_roi%d_%d_nsample_%d zeta lP_record xi option',idx_roi,ith_run,Nsubsample));

saveid = 1;
