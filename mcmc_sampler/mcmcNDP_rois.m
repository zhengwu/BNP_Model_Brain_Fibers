%for nDP process
function [zeta, lP_record, xi, saveid] = mcmcNDP_rois(ith_run,idx_roi,In_COMP2,COMP1,ALPHA_BETA_PRIOR)
    % input: ith_run - index for this run   
    %        idx_roi - index for rois, taking the value from 1 to 45
    %        (according to our data)
    
    %        In_COMP2 - indicator for including the translation component
    %        COMP1 - name of component 1, can have two values,'fpca' or 'trans'
    %        ALPHA_BETA_PRIOR - indicator to set priors for alpha and beta
    
    addpath(genpath('./distributions'));
    addpath(genpath('./utilities'));
    addpath(genpath('./curve_processing_functions'));
    
    %parameters
    Npoint = 100;
    Kcluster = 1;
    
    %parameter for subsampling
    SubSamp =1;
    Nsubsample = 400;

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
   
    %%%%%%%%%%%%%%%%%%% process rotation %%%%%%%%%%%%%%%%%%%%
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
        inputdata{idx}.Rots = mapedR;
        clear mapedR;
    end
    
    option.INDICATOR_COMP3=0;
    option.INDICATOR_COMP2=In_COMP2;
    option.COMP1 = COMP1;
    option.NIW = 1;
    option.GRAPHICS=0;
    option.ALPHA_BETA_PRIOR = ALPHA_BETA_PRIOR;
    option.FIRST_LS = 1; %first label-swithing;
    option.SECOND_LS = 1; %second label-swithing;

    %[zeta,lP_record,xi] = samplerNDP(inputdata,9,15,9,15,5000,option);
    [zeta,lP_record,xi] = samplerNDP(inputdata,9,15,9,15,5000,option);
    
    eval(sprintf('save norot_mcmcresult_roi%d_%d_nsample_%d zeta lP_record xi option',idx_roi,ith_run,Nsubsample));

saveid = 1;
