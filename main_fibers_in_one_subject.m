%for modeling fiber curves between two regions of interet in one subject

clear all; 
close all;

% %%%%%%%%%%%%%%%%%%% add necessary paths %%%%%%%%%%%%%%%%%%%%
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
     addpath(genpath('./distributions'))
end

% %%%%%%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% dataset 1 %%%%%%%%%%
roia = 16;
roib = 55;
isub = 1;
iscan = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%% dataset 2 %%%%%%%%%%
% roia = 16;
% roib = 56;
% isub = 1;
% iscan = 1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kcluster = 1;
Npoint = 100;
pre = 'desikan';
eval(sprintf('load R_test_retest-atlas_data_%s_K%d_roi%d-%d_%d',pre,Kcluster,roia,roib,Npoint));
pca_basis = pca_result.pca_basis;
load(sprintf('Rtest-retest_A%d_S%d_%s_roi%d-%d_streamlines_fPCAresults_%d',isub,iscan,pre,roia,roib,Npoint));

Nfiber = length(labels);
newfibers = normal_fibers;

%process rotation;
Omean = Mean_Rotation(Rots);
for i=1:size(Rots,3)
    newRots(:,:,i) = Omean'*Rots(:,:,i);
end
Rots = newRots;

N = Nfiber;
fpca_coeff = [];
for i=1:Nfiber
    id = i;
    tmp_coeff = PCA_coeffs{id};
    coeff_vect = [tmp_coeff{1}(1),tmp_coeff{2}(1),tmp_coeff{3}(1)]';
    fpca_coeff(:,i) = coeff_vect;
    new_PCA_coeffs{i} = tmp_coeff;
end

for i=1:N
    R = Rots(:,:,i);
    trR = trace(R);
    theta(i) = acos((trR-1)/2);
    tempM = (R - R')*theta(i)/(2*sin(theta(i)));
    mapedR(:,i) = [-R(1,2),R(1,3),R(3,2)]';
end

%reconstruction
final_reconstl = roifibers_reconstruction_adaptive(labels,Trans,Rots,new_PCA_coeffs,means_sl,pca_basis);
clear new_PCA_coeffs;

%reconstruction quality check
dt = 1/(Npoint-1);
point_off = zeros(1,N);
for j=1:N
    point_off(j)=sqrt(sum((sum((final_reconstl(:,:,j)-normal_fibers(:,:,j)).^2,1)*dt)));
end

for j=1:N
    for i=1:size(final_reconstl,2)
        C(j,i,1) = norm(final_reconstl(:,i,j)-newfibers(:,i,j));
        C(j,i,2) = norm(final_reconstl(:,i,j)-newfibers(:,i,j));
    end
end

figure(10*isub+iscan);clf; hold on;
for j=1:N
    mesh([squeeze(final_reconstl(1,:,j))',squeeze(final_reconstl(1,:,j))'],[squeeze(final_reconstl(2,:,j))',squeeze(final_reconstl(2,:,j))'],[squeeze(final_reconstl(3,:,j))',squeeze(final_reconstl(3,:,j))'],squeeze(C(j,:,:)));
end
title('Reconstructed Streamlines');
axis off;
set(gca,'fontsize',22);
view([0,76]);
  
figure(100*isub+iscan);clf; hold on;
for j=1:N
    plot3(newfibers(1,:,j),newfibers(2,:,j),newfibers(3,:,j),'linewidth',3);
end
title('Original Streamlines');
axis off;
set(gca,'fontsize',22);
view([0,76]);


%plot different components
figure(8);clf;
color = {'ko', 'ro', 'go', 'bo', 'mo', 'co'};
hold on;
scatter3(fpca_coeff(1,:),fpca_coeff(2,:),fpca_coeff(3,:),color{min(isub)});
set(gca,'fontsize',22)
title('Shape component');

figure(9);clf;
color = {'ko', 'ro', 'go', 'bo', 'mo', 'co'};
hold on;
scatter3(Trans(1,:),Trans(2,:),Trans(3,:),color{isub});
set(gca,'fontsize',22)
title('Trans component');


figure(10);clf;
color = {'ko', 'ro', 'go', 'bo', 'mo', 'co'};
hold on;
scatter3(mapedR(1,:),mapedR(2,:),mapedR(3,:),color{isub});
set(gca,'fontsize',22)
title('Rots component');


%%%%%%%%%%%%%%%%%%% get the data ready for modeling %%%%%%%%%%%%%%%%%%%%
%%% here we perform normalization 

comp1 = fpca_coeff;
comp2 = Trans;
comp3 = mapedR;

mean_comp1 = mean(comp1,2);
mean_comp2 = mean(comp2,2);
mean_comp3 = mean(comp3,2);

%remove the mean
comp1_c = comp1 - mean_comp1*ones(1, size(comp1,2));
comp2_c = comp2 - mean_comp2*ones(1, size(comp2,2));
comp3_c = comp3 - mean_comp3*ones(1, size(comp3,2));

%rescale the data;
for i=1:3
    varc1(i) = var(comp1(i,:));
    varc2(i) = var(comp2(i,:));
    varc3(i) = var(comp3(i,:));
    
    comp1_c(i,:) = comp1_c(i,:)/sqrt(varc1(i));
    comp2_c(i,:) = comp2_c(i,:)/sqrt(varc2(i));
    comp3_c(i,:) = comp3_c(i,:)/sqrt(varc3(i));
end

inputdatan.fpca_coeff = comp1_c;
inputdatan.Trans = comp2_c;
inputdatan.mapedR = comp3_c;
inputdatan.Rots = Rots;
inputdatan.Reconfiber = final_reconstl;

%keyboard;
close all; %close all figures



Nsample = 10000;

% seed for reproducibility

%d=rng;
%save random_seed d;
load('random_seed.mat');
rng(d);

tic
% method 1 
% sample based on CRP - faster
%[class_id, K_record, components, alpha_record] = sampler_crp(inputdatan, Nsample);

% method 2 - reported in the paper
% sample based on Dirichlet
[class_id, K_record, components, alpha_record] = sampler_dir(inputdatan,10,Nsample);

%change to use different components for modeling
% %%%%%%%%%%%%%%%%%%% get posterior samples %%%%%%%%%%%%%%%%%%%%
% to use different combination of components, please edit parameters inside
% the function "sampler_crp" or "sampler_dir". For example, 
% (1) if use only the shape component to perform the modeling, set
%      INDICATOR_COMP1 = 1; 
%      INDICATOR_COMP2 = 0; 
%      INDICATOR_COMP3 = 0; 
%      component1 = fpca_coeff;
%      component2 = trans;
%      component3 = mappedrots;

% (2) if use only the translation component to perform the modeling, set
%      INDICATOR_COMP1 = 1;
%      INDICATOR_COMP2 = 0; 
%      INDICATOR_COMP3 = 0;
%      component2 = fpca_coeff;
%      component1 = trans;
%      component3 = mappedrots;


% (3) if we want to all components to perform the modeling, set
%      INDICATOR_COMP1 = 1; 
%      INDICATOR_COMP2 = 1; 
%      INDICATOR_COMP3 = 1;
%      component1 = fpca_coeff;
%      component2 = trans;
%      component3 = mappedrots;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


toc


% %%%%%%%%%%%%%%%%%%% show posterior results %%%%%%%%%%%%%%%%%%%%

record_label = class_id';
burnin = 1000;
N = size(record_label,2);
%plot the original data;
finalB = zeros(N,N);
for i=burnin:Nsample
   tempB = calculate_B(squeeze(record_label(i,:)),N);
   finalB = finalB + tempB;
   psoteriorK(i-burnin+1) = length(unique(record_label(i,:)));
end

%show the posterior distribution of K
figure(98);clf;
histogram(psoteriorK);
set(gca,'fontsize',22);

%compare to the "true" label;
ind_modK = mode(psoteriorK);
[c,cn,thrd] = BtoCluster_zw(finalB,N,mode(psoteriorK));
display('The final clustering result is:');
display(c);
eval(sprintf('load roi%d_%d_sub%d_scan%d_truelable.mat',roia,roib,isub,iscan));
[AR,RI,MI,HI] = RandIndex(c,ture_c)


%show pairwise probability matrix
D_final = [];
M = length(unique(ture_c));
cct(1) = 1;
new_order = [];
for i=1:M
    class{i} = find(ture_c==i);
    cct(i+1) = cct(i) +  length(class{i})-1;
    new_order = [new_order; class{i}];
end

for i=1:N
    for j=i+1:N
        D_final(i,j)= finalB(new_order(i),new_order(j));
        D_final(j,i)=D_final(i,j);
    end
end
figure(99);clf; imagesc(D_final/max(max(D_final)));
set(gca,'fontsize',20);
colorbar;


%plots some component with clustering information
% comp =  comp3_c;
% 
% color = {'ro', 'go', 'bo', 'mo', 'co', 'yo'};
% figure(101);clf;
% for i=1:min(ind_modK,7)
%     cluster_idx = find(c==i);
%     for j = 1:length(cluster_idx)
%         hold on;
%         scatter3(comp(1,cluster_idx(j)),comp(2,cluster_idx(j)),comp(3,cluster_idx(j)),color{i});
%     end
% end
% set(gca,'fontsize',22);
% title('Final clustering results')

%plot the fibers
color = { 'r', 'g', 'b', 'm', 'c','y','k'};
figure(102);clf;
for i=1:min(ind_modK,7)
    cluster_idx = find(c==i);
    for j = 1:length(cluster_idx)
        hold on;
        plot3(newfibers(1,:,cluster_idx(j)),newfibers(2,:,cluster_idx(j)),newfibers(3,:,cluster_idx(j)),color{i});
    end 
end
%view([-5,18]);
view([-115,10]);
set(gca,'fontsize',22);
title('Final clustering results')