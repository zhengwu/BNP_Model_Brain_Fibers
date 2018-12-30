function [zeta,lP_record,xi] = samplerNDP_rotation(input_data,K,L,intK,intL,num_sweeps,option)
% main function for sampling posteriors from mixture of product kernels

INDICATOR_COMP3 = option.INDICATOR_COMP3;
INDICATOR_COMP2 = option.INDICATOR_COMP2;
NIW =     option.NIW;
GRAPHICS=    option.GRAPHICS;
ALPHA_BETA_PRIOR=   option.ALPHA_BETA_PRIOR;
FIRST_LS = option.FIRST_LS;  % first label-swithing;
SECOND_LS=    option.SECOND_LS; %second label-swithing;

    
    
%get the data from the input data;
NSUB = length(input_data);

%align the data by remove the individual mean
for isub = 1:NSUB
    tran_submean = mean(input_data{isub}.Trans,2);
    input_data{isub}.Trans = input_data{isub}.Trans-tran_submean*ones(1,size(input_data{isub}.Trans,2));
end

% center the data;
comp1all = [];
comp2all = [];
allidx = 0;
for isub = 1:NSUB
     comp1all = [comp1all,input_data{isub}.fpca_coeff];
     comp2all = [comp2all,input_data{isub}.Trans];
     allidx = allidx + size(input_data{isub}.fpca_coeff,2);
end

mean_comp1 = mean(comp1all,2);
mean_comp2 = mean(comp2all,2);

%remove the global mean
for isub = 1:NSUB
     input_data{isub}.fpca_coeff = input_data{isub}.fpca_coeff-mean_comp1*ones(1,size(input_data{isub}.fpca_coeff,2));
     input_data{isub}.Trans = input_data{isub}.Trans-mean_comp2*ones(1,size(input_data{isub}.Trans,2));
end

% rescale the data;
for i=1:3
   varc1(i) = var(comp1all(i,:));
   varc2(i) = var(comp2all(i,:));
end

for isub = 1:NSUB
for i=1:3
     input_data{isub}.fpca_coeff(i,:) = input_data{isub}.fpca_coeff(i,:)./sqrt(varc1(i));
     input_data{isub}.Trans(i,:) = input_data{isub}.Trans(i,:)./sqrt(varc2(i));
end
end

%plot the data and check;
color = {'ko', 'ro', 'go', 'bo', 'mo', 'co','k*','r*','b*','m*','c*','k.','r.','b.'};
figure(10);clf;
figure(9);clf;
for isub = 1:NSUB
    idsub = ceil(isub/3);
    figure(1);
    %plot the shape part;
    figure(10);
    hold on;
    scatter3(input_data{isub}.fpca_coeff(1,:),input_data{isub}.fpca_coeff(2,:),input_data{isub}.fpca_coeff(3,:),color{idsub});
    
    
    figure(9);
    hold on;
    scatter3(input_data{isub}.Trans(1,:),input_data{isub}.Trans(2,:),input_data{isub}.Trans(3,:),color{idsub});
end


for isub = 1:NSUB
    %pca coefficent
    tmp_fpca_coeff = input_data{isub}.fpca_coeff;
    %tmp_fpca_coeff = tmp_fpca_coeff - mean_comp1*ones(1,size(tmp_fpca_coeff,2));
    fpca_coeff{isub} = tmp_fpca_coeff;
    
    %translation
    tmp_Trans = input_data{isub}.Trans;
    %tmp_Trans = tmp_Trans - mean_comp2*ones(1,size(tmp_Trans,2));
    trans{isub} = tmp_Trans;
    
    %original rotations
    tmp_Rots = input_data{isub}.Rots;
    rots{isub} = input_data{isub}.Rots;
    
    %mapped rotations
    tmp_mappedrots = input_data{isub}.mapedR;
   % tmp_mappedrots = tmp_mappedrots - mean(tmp_mappedrots,2)*ones(1,size(tmp_mappedrots,2));
    mappedrots{isub} = tmp_mappedrots;
    
    % number of curves in each connection
    ObjN(isub) = size(input_data{isub}.Trans,2);
end


if(strcmp(option.COMP1, 'fpca'))
    component1 = fpca_coeff;
    component2 = trans;
else
    component2 = fpca_coeff;
    component1 = trans;
end
component3 = mappedrots;


if(nargin < 2)
    num_sweeps = 1000;
end


%%%%%%%%%%%%%%%% initialize parameters%%%%%%%%%%%%%%%%
% initialize center indicators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for nDP
K = K;
L = L;

%for nDP
%initial value for alpha and beta
alpha = 1;
beta = 1;

%prior on alpha and beta;
a_alpha = 3;
b_alpha = 3;

a_beta = 3;
b_beta = 3;


% set alpha's gamma prior parameters
a_0 = 1;
b_0 = 1;


D1 = size(component1{1},1);
D2 = size(component2{1},1);
D3 = size(component3{1},1);

%%%%%%%%%%%%%%%% set the hyper parameters%%%%%%%%%%%%%%%%
% set normal inverse wishart hyper parameters and gamma prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set normal inverse wishart hyper parameters for comp1
mu_0_comp1 = zeros(size(component1{1}(:,1)));
k_0_comp1 = 1;
lambda_0_comp1 = eye(3,3);%cov(component1{1}(:,1:10)');
v_0_comp1 = 5;

% set normal inverse wishart hyper parameters for comp2
mu_0_comp2 = zeros(size(component2{1}(:,1)));
k_0_comp2 = 1;
lambda_0_comp2 = eye(3,3);%cov(component2{1}(:,1:10)');
v_0_comp2 = 6;

% set hyper parameters for comp3
mu_0_comp3 = zeros(size(component3{1}(:,1)));
k_0_comp3 = 1;
v_0_comp3 = size(component3{1}(:,1),1);
lambda_0_comp3 = eye(3,3);%cov(component3{1}(:,1:50)');

%for center indicators
zeta = zeros(num_sweeps,NSUB);
pi_ci = zeros(num_sweeps,K);
for isub = 1:NSUB
    zeta(1,isub) = 1+mod(isub,intK);
end;

for i = 1:K
    pi_ci(1,i) = 1/K;
end;

%for group indicators
xi{num_sweeps,NSUB} = {};
for j=1:NSUB
    nobj = ObjN(j);
    tmpidx = 1+mod(1:nobj,intL);  
    xi{1,j} = tmpidx;
end

%group weights
w = zeros(num_sweeps,K,L);
for k=1:K
    for l=1:L
        w(1,k,l) = 1/L;
    end;
end;


%log probability record
lP_record = zeros(num_sweeps,1);


%for each component, prepare the space for recording their distribution
%parameters

for k=1:K
    for l=1:L
        %for component1
        mean_record_comp1{k,l} = cell(num_sweeps,1);
        covariance_record_comp1{k,l} = cell(num_sweeps,1);
        inv_covariance_record_comp1{k,l} = cell(num_sweeps,1);
        
        %for component2
        mean_record_comp2{k,l} = cell(num_sweeps,1);
        covariance_record_comp2{k,l} = cell(num_sweeps,1);
        inv_covariance_record_comp2{k,l} = cell(num_sweeps,1);
        
        %for component3
        mean_record_comp3{K,l} = cell(num_sweeps,1);
        covariance_record_comp3{k,l} = cell(num_sweeps,1);
        inv_covariance_record_comp3{k,l} = cell(num_sweeps,1);
    end
end


%for each component, initialize the parameters
for k=1:K
    subK_ind = find(zeta(1,:)==k);
    if(length(subK_ind)==0)
        %assign priors
        for l=1:L
            %for component 1
            new_covariance_comp1 = iwishrnd(inv(lambda_0_comp1), v_0_comp1+3);
            mean_record_comp1{k,l}{1} = mvnrnd(mu_0_comp1',new_covariance_comp1'/k_0_comp1)';%component1{1+mod(k,NSUB)}(:,l);
            covariance_record_comp1{k,l}{1} = new_covariance_comp1;
            inv_covariance_record_comp1{k,l}{1} = inv(covariance_record_comp1{k,l}{1});
            
            %for component 2
            new_covariance_comp2 = iwishrnd(inv(lambda_0_comp2), v_0_comp2+3);
            mean_record_comp2{k,l}{1} = mvnrnd(mu_0_comp2',new_covariance_comp2'/k_0_comp2)';%component1{1+mod(k,NSUB)}(:,l);
            covariance_record_comp2{k,l}{1} = new_covariance_comp2;
            inv_covariance_record_comp2{k,l}{1} = inv(covariance_record_comp2{k,l}{1});
            
        end
    else
     for l=1:L
         subK_data_likelihood_lp = zeros(size(subK_ind));
         
         kl_comp1 = zeros(D1,1);
         kl_comp2 = zeros(D2,1);
         kl_comp3 = zeros(D3,1);
         tmp_idx = 1;
         for i=1:length(subK_ind)
             isub = subK_ind(i);
             obj_assign = xi{1,isub};
             %get the data
             curr_comp1 = component1{isub};
             curr_comp2 = component2{isub};
             curr_comp3 = component3{isub};

             clear obj_l;
             obj_l = find(obj_assign==l);
             
             for iobj=1:length(obj_l)
                 objidx = obj_l(iobj);
                 kl_comp1(:,tmp_idx) = curr_comp1(:,objidx);
                 kl_comp2(:,tmp_idx) = curr_comp2(:,objidx);
                 kl_comp3(:,tmp_idx) = curr_comp3(:,objidx);
                 tmp_idx = tmp_idx + 1;
             end;
             
         end
         mean_record_comp1{k,l}{1} = mean(kl_comp1,2);
         covariance_record_comp1{k,l}{1} = cov(kl_comp1') + 0.1*eye(3,3);
         inv_covariance_record_comp1{k,l}{1} = inv(covariance_record_comp1{k,l}{1});
         
         mean_record_comp2{k,l}{1} = mean(kl_comp2,2);
         covariance_record_comp2{k,l}{1} = cov(kl_comp2') + 0.1*eye(3,3);
         inv_covariance_record_comp2{k,l}{1} = inv(covariance_record_comp2{k,l}{1});
         
    end
end
end


for k=1:K
    for l=1:L
        %for component 3
        new_covariance_comp3 = iwishrnd(lambda_0_comp3, v_0_comp3);
        mean_record_comp3{k,l}{1} = mvnrnd(mu_0_comp3',new_covariance_comp3'/k_0_comp3)';
        covariance_record_comp3{k,l}{1} = new_covariance_comp3;%eye(size(component3{1},1));
        inv_covariance_record_comp3{k,l}{1} = eye(size(component3{1},1));
    end
end


% run the Gibbs sampler
for sweep = 2:num_sweeps 
    
    disp(['Sweep ' num2str(sweep) '/' num2str(num_sweeps)])
    
    
    for k=1:K
        for l=1:L
            mean_record_comp1{k,l}{sweep} = mean_record_comp1{k,l}{sweep-1};
            covariance_record_comp1{k,l}{sweep} = covariance_record_comp1{k,l}{sweep-1};
            inv_covariance_record_comp1{k,l}{sweep} = inv_covariance_record_comp1{k,l}{sweep-1};
            
            if(INDICATOR_COMP2)
                mean_record_comp2{k,l}{sweep} = mean_record_comp2{k,l}{sweep-1};
                covariance_record_comp2{k,l}{sweep} = covariance_record_comp2{k,l}{sweep-1};
                inv_covariance_record_comp2{k,l}{sweep} = inv_covariance_record_comp2{k,l}{sweep-1};
            end
            
            if(INDICATOR_COMP3)
                mean_record_comp3{k,l}{sweep} = mean_record_comp3{k,l}{sweep-1};
                covariance_record_comp3{k,l}{sweep} = covariance_record_comp3{k,l}{sweep-1};
                inv_covariance_record_comp3{k,l}{sweep} = inv_covariance_record_comp3{k,l}{sweep-1};
            end
    
        end
    end
    
   
   %center indicator
   zeta(sweep,:) = zeta(sweep-1,:);
   pi_ci(sweep,:) = pi_ci(sweep-1,:);
   
   %grop indicator
   for isub = 1:NSUB
       xi{sweep,isub} = xi{sweep-1,isub};
   end
   w(sweep,:,:) = w(sweep-1,:,:);
   
   
   for isub = 1:NSUB
       %for each subject, reassign a cluster to it
       
       %get the data;
       curr_comp1 = component1{isub};
       curr_comp2 = component2{isub};
       curr_comp3 = component3{isub};
       
       %calculate the probability assigning each subject to a class
       data_likelihood_lp = zeros(K,1);
       for k=1:K
           %prior
           tmp_pi = pi_ci(sweep,k);
           atoms_likelihood = zeros(ObjN(isub),L);
           for Ij = 1:ObjN(isub) %each element in the class
               for l=1:L
                   %likelihood
                   l1 = curr_comp1(:,Ij) - mean_record_comp1{k,l}{sweep};
                   if(INDICATOR_COMP2)
                       l2 = curr_comp2(:,Ij) - mean_record_comp2{k,l}{sweep};
                   end
                   if(INDICATOR_COMP3)
                       l3 = curr_comp3(:,Ij) - mean_record_comp3{k,l}{sweep};
                   end
                   
                   %likelihood for each component;
                   lp1 = -l1'*inv_covariance_record_comp1{k,l}{sweep}*l1/2-.5*log(det(covariance_record_comp1{k,l}{sweep}))-(D1/2)*log(2*pi);
                   if(INDICATOR_COMP2)
                       lp2 = -l2'*inv_covariance_record_comp2{k,l}{sweep}*l2/2-.5*log(det(covariance_record_comp2{k,l}{sweep}))-(D2/2)*log(2*pi);
                   end
                   if(INDICATOR_COMP3)
                    %lp3 = LogLikelihoodDataTerm(mean_record_comp3{k,l}{sweep},kappa_comp3{sweep}(l),l3);
                     lp3 = -l3'*inv_covariance_record_comp3{k,l}{sweep}*l3/2-.5*log(det(covariance_record_comp3{k,l}{sweep}))-(D3/2)*log(2*pi);
                   end
            
                   if(INDICATOR_COMP3&&INDICATOR_COMP2)
                       data_likelihood_kl_lp = lp1+lp2+lp3;
                   elseif(INDICATOR_COMP2)
                       data_likelihood_kl_lp = lp1+lp2;
                   else
                       data_likelihood_kl_lp = lp1;
                   end
                   
                   %prior
                   tmp_w_kl = w(sweep,k,l); %weight for (k,l);
                   
                   atoms_likelihood(Ij,l) = log(tmp_w_kl)+(data_likelihood_kl_lp);
               end
           end
           if(sweep>2)
               maxe = 0;
               %maxe = max(max(atoms_likelihood));
               atoms_likelihood = atoms_likelihood - maxe;
               tmp_likelihood = sum(exp(atoms_likelihood),2);
               %approximate
               %tmp_likelihood(find(tmp_likelihood==0)) = min(tmp_likelihood(find(tmp_likelihood>0)));
           else
                maxe = 0;
%                tmp_likelihood = sum(exp(atoms_likelihood),2);
                %maxe = max(max(atoms_likelihood));
                atoms_likelihood = atoms_likelihood - maxe;
                tmp_likelihood = sum(exp(atoms_likelihood),2);
               %approximate
               %tmp_likelihood(find(tmp_likelihood==0)) = min(tmp_likelihood(find(tmp_likelihood>0)));
           end
           
           data_likelihood_lp(k) = log(tmp_pi)+sum((log(tmp_likelihood)+maxe));
           test(k) = sum(log(tmp_likelihood));
       end
       
       infidx = isinf(data_likelihood_lp);
       data_likelihood_lp(infidx)=-10^10;
       data_likelihood_lp = data_likelihood_lp - max(data_likelihood_lp); 
       data_likelihood = exp(data_likelihood_lp);
       data_likelihood(infidx) = 0;
       
       if(sum(data_likelihood)==0)
           data_likelihood = rand(size(data_likelihood));
           display('randomly choose data_likelihood')
       end
       data_likelihood = data_likelihood/sum(data_likelihood);
       
       cdf = cumsum(data_likelihood);
       rn = rand;
       
       zeta(sweep,isub) = min(find((cdf>rn)==1));
   end

   %%%%%%%%%%%%% step 3, sample pi_k%%%%%%%%%%%%%
   tmp_zeta = zeta(sweep,:);
   for k=1:K
       m(k) = length(find(tmp_zeta==k));
   end
   
   u = zeros(K,1);
   for k=1:K-1
       u(k)= betarnd(1+m(k), alpha + sum(m(k+1:end)));
   end;
   u(K) = 1;
   
   for k = 1:K
       pi_ci(sweep,k) = u(k)*prod(1-u(1:k-1));
   end;
 
 
%    %%%%%%%%% propose random changes%%%%%%%
   %propose random changes every KK iter
   
   if(sweep>2000)
       KK = 20;
   else
       KK = 2;
   end
   
   if(mod(sweep,KK)==0&&FIRST_LS==1)
       % first type of label-switching
       if(NSUB>1)
           all_labels = unique(zeta(sweep,:));
           perm_label = randperm(length(all_labels));
           if(length(perm_label)>1)
               idx1 = find(zeta(sweep,:) == all_labels(perm_label(1)));
               idx2 = find(zeta(sweep,:) == all_labels(perm_label(2)));
               label1 = all_labels(perm_label(1));
               label2 = all_labels(perm_label(2));
               
               %decide whether accept this label
               p1 = pi_ci(sweep,label1);
               p2 = pi_ci(sweep,label2);
               m1 = length(idx1);
               m2 = length(idx2);
               r = min([1,(p1/p2)^(m2-m1)]);
               rand_u = rand;
               if(rand_u<r) %accept this swith
                   zeta(sweep,idx1(randi([1,length(idx1)],1))) = label2;
                   zeta(sweep,idx2(randi([1,length(idx2)],1))) = label1;
                   
                   display('label swithed')
               end
           end
       end
       
       % second type of label-switching
       if(NSUB>1&&SECOND_LS==1)
           selected_k = min(K-1,ceil(rand()*K));
           %acceptance rate
           a_rate = min(1, (1-u(selected_k+1))^zeta(sweep,selected_k)/(1-u(selected_k))^zeta(sweep,selected_k+1) );
           rand_u = rand;
           
           %accept changes
           if(rand_u<a_rate)
               tmp = u(selected_k);
               u(selected_k) = u(selected_k+1);
               u(selected_k+1) = tmp;
           end
           %update pi_ci
           for k = 1:K
               pi_ci(sweep,k) = u(k)*prod(1-u(1:k-1));
           end
       end
       
   end
   

   %%%%%%%%%%%%% step 2, sample group indicator%%%%%%%%%%%%%
   for isub = 1:NSUB
       %get the data
       curr_comp1 = component1{isub};
       curr_comp2 = component2{isub};
       curr_comp3 = component3{isub};
       
       for iobj = 1:ObjN(isub)
           itom_likelihood_lp = zeros(L,1);
           for l = 1:L
               
               %data likelihood
               l1 = curr_comp1(:,iobj) - mean_record_comp1{zeta(sweep,isub),l}{sweep};
               if(INDICATOR_COMP2)
                   l2 = curr_comp2(:,iobj) - mean_record_comp2{zeta(sweep,isub),l}{sweep};
               end
               if(INDICATOR_COMP3)
                   l3 = curr_comp3(:,iobj) - mean_record_comp3{zeta(sweep,isub),l}{sweep};
               end
               
               %likelihood for each component;
               lp1 = -l1'*inv_covariance_record_comp1{zeta(sweep,isub),l}{sweep}*l1/2-.5*log(det(covariance_record_comp1{zeta(sweep,isub),l}{sweep}))-(D1/2)*log(2*pi);
               
               if(INDICATOR_COMP2)
                   lp2 = -l2'*inv_covariance_record_comp2{zeta(sweep,isub),l}{sweep}*l2/2-.5*log(det(covariance_record_comp2{zeta(sweep,isub),l}{sweep}))-(D2/2)*log(2*pi);
               end
               
               if(INDICATOR_COMP3)
                   %lp3 = LogLikelihoodDataTerm(mean_record_comp3{k,l}{sweep},kappa_comp3{sweep}(l),l3);
                   lp3 = -l3'*inv_covariance_record_comp3{zeta(sweep,isub),l}{sweep}*l3/2-.5*log(det(covariance_record_comp3{zeta(sweep,isub),l}{sweep}))-(D3/2)*log(2*pi);
               end
               
               if(INDICATOR_COMP3&&INDICATOR_COMP2)
                   data_likelihood_kl_lp = lp1+lp2+lp3;
               elseif(INDICATOR_COMP2)
                   data_likelihood_kl_lp = lp1+lp2;
               else
                   data_likelihood_kl_lp = lp1;
               end
               
               %prior
               tmp_w_kl = w(sweep,zeta(sweep,isub),l); %weight for (k,l);
               
               itom_likelihood_lp(l) = log(tmp_w_kl)+(data_likelihood_kl_lp);
           end
           % sample the group indicators
           
           infidx = isinf(itom_likelihood_lp);
           itom_likelihood_lp(infidx)=-10^10;
           itom_likelihood_lp = itom_likelihood_lp - max(itom_likelihood_lp);
           itom_likelihood = exp(itom_likelihood_lp);
           itom_likelihood(infidx) = 0;
       
           if(sum(itom_likelihood)==0)
               itom_likelihood(end) = eps;
               display('random sample itom_likelihood')
           end
           itom_likelihood = itom_likelihood/sum(itom_likelihood);
                  
           cdf = cumsum(itom_likelihood);
           rn = rand;
           xi{sweep,isub}(iobj) =  min(find((cdf>rn)==1));
       end
   end;
   
    
   
 %%%%%%%%%%%%% step 4, sample w_lk%%%%%%%%%%%%%
 tmp_zeta = zeta(sweep,:);
 n = zeros(K,L);
 for k=1:K
     subK_ind = find(tmp_zeta==k);
     for i=1:length(subK_ind)
         isub = subK_ind(i);
         for iobj = 1:ObjN(isub)
             objidx = xi{sweep,isub}(iobj);
             n(k,objidx) = n(k,objidx) + 1;
         end
     end;
 end
 figure(1);
 imagesc(n);
 pause(0.01);
 
 v = zeros(K,L);
 for k=1:K
     for l=1:L-1
          v(k,l) = betarnd(1+n(k,l), beta+sum(n(k,l+1:end),2));
%           v(k,l) = betarnd(1+n(k,l), beta+sum(n(k,l+1:end),2));
     end
     v(k,L) = 1;
 end
 
 for k=1:K
     for l=1:L
        w(sweep,k,l) = v(k,l)*prod(1 - v(k,1:l-1)); 
     end
 end

 
 %%%%%%%%%%%%% step 5, sample \theta_lk%%%%%%%%%%%%%
 tic
  for k=1:K
     subK_ind = find(tmp_zeta==k);
     for l=1:L
         subK_data_likelihood_lp = zeros(size(subK_ind));
         
         kl_comp1 = zeros(D1,1);
         kl_comp2 = zeros(D2,1);
         kl_comp3 = zeros(D3,1);
         tmp_idx = 1;
         for i=1:length(subK_ind)
             isub = subK_ind(i);
             obj_assign = xi{sweep,isub};
             %get the data
             curr_comp1 = component1{isub};
             curr_comp2 = component2{isub};
             curr_comp3 = component3{isub};

             clear obj_l;
             obj_l = find(obj_assign==l);
             
             for iobj=1:length(obj_l)
                 objidx = obj_l(iobj);
                 kl_comp1(:,tmp_idx) = curr_comp1(:,objidx);
                 kl_comp2(:,tmp_idx) = curr_comp2(:,objidx);
                 kl_comp3(:,tmp_idx) = curr_comp3(:,objidx);
                 tmp_idx = tmp_idx + 1;
             end;
             
         end
         
         %posterior for theta_lk;
         
         N_this_class = tmp_idx - 1;
         
         if(N_this_class==0) % if there is no elements in a class, we sample from the prior distribution
             new_covariance_comp1 = iwishrnd(lambda_0_comp1, v_0_comp1);
             new_mean_comp1 = mvnrnd(mu_0_comp1',new_covariance_comp1'/k_0_comp1)';
             
             new_covariance_comp2 = iwishrnd(lambda_0_comp2, v_0_comp2);
             new_mean_comp2 = mvnrnd(mu_0_comp2',new_covariance_comp2'/k_0_comp2)';

             new_covariance_comp3 = iwishrnd(lambda_0_comp3, v_0_comp3);
             new_mean_comp3 = mvnrnd(mu_0_comp3',new_covariance_comp3'/k_0_comp3)';
             
             
         else
             if(N_this_class>3)
                 S1 = cov(kl_comp1')'*(N_this_class-1);
                 y_bar1 = mean(kl_comp1,2);
                 mu_n_comp1 = (N_this_class*y_bar1+...
                     k_0_comp1*mu_0_comp1)/(k_0_comp1+N_this_class);
             else
                 y_bar1 = mean(kl_comp1,2);
                 S1 = (kl_comp1 - repmat(y_bar1,1,size(kl_comp1,2)))*(kl_comp1 - repmat(y_bar1,1,size(kl_comp1,2)))';
                 mu_n_comp1 = (N_this_class*y_bar1+...
                     k_0_comp1*mu_0_comp1)/(k_0_comp1+N_this_class);
             end
             
             if(INDICATOR_COMP2)
                 
                 if(N_this_class>3)
                     S2 = cov(kl_comp2')'*(N_this_class-1);
                     y_bar2 = mean(kl_comp2,2);
                     mu_n_comp2 = (N_this_class*y_bar2+...
                         k_0_comp2*mu_0_comp2)/(k_0_comp2+N_this_class);
                 else
                     y_bar2 = mean(kl_comp2,2);
                     S2 = (kl_comp2 - repmat(y_bar2,1,size(kl_comp2,2)))*(kl_comp2 - repmat(y_bar2,1,size(kl_comp2,2)))';
                     mu_n_comp2 = (N_this_class*y_bar2+...
                         k_0_comp2*mu_0_comp2)/(k_0_comp2+N_this_class);
                 end
             
             end
             
             if(INDICATOR_COMP3)
                 S3 = cov(kl_comp3')'*(N_this_class-1);
                 y_bar3 = mean(kl_comp3,2);
                 mu_n_comp3 = (N_this_class*y_bar3+...
                     k_0_comp3*mu_0_comp3)/(k_0_comp3+N_this_class);
             end
             
        
             if(NIW)
                 
                 lambda_n_comp1 = lambda_0_comp1 + ...
                     (k_0_comp1*N_this_class)/(k_0_comp1+N_this_class)*...
                     ((y_bar1-mu_0_comp1)*(y_bar1-mu_0_comp1)')+S1;
                 
                 k_n_comp1 = k_0_comp1 + N_this_class;
                 v_n_comp1 = v_0_comp1 + N_this_class;
                 % sample the new covariance and mean from the joint posterior
                 new_covariance_comp1 = iwishrnd(lambda_n_comp1, v_n_comp1);
                 new_mean_comp1 = mvnrnd(mu_n_comp1',new_covariance_comp1'/k_n_comp1)';
                 
                 if(INDICATOR_COMP2)
                     lambda_n_comp2 = lambda_0_comp2 + ...
                         (k_0_comp2*N_this_class)/(k_0_comp2+N_this_class)*...
                         ((y_bar2-mu_0_comp2)*(y_bar2-mu_0_comp2)')+S2;
                     
                     k_n_comp2 = k_0_comp2 + N_this_class;
                     v_n_comp2 = v_0_comp2 + N_this_class;
                     % sample the new covariance and mean from the joint posterior
                     new_covariance_comp2 = iwishrnd(lambda_n_comp2, v_n_comp2);
                     new_mean_comp2 = mvnrnd(mu_n_comp2',new_covariance_comp2'/k_n_comp2)';
                 end
                 
                 if(INDICATOR_COMP3)
                     lambda_n_comp3 = lambda_0_comp3 + ...
                         (k_0_comp3*N_this_class)/(k_0_comp3+N_this_class)*...
                         ((y_bar3-mu_0_comp3)*(y_bar3-mu_0_comp3)')+S3;
                     
                     k_n_comp3 = k_0_comp3 + N_this_class;
                     v_n_comp3 = v_0_comp3 + N_this_class;
                     % sample the new covariance and mean from the joint posterior
                     new_covariance_comp3 = iwishrnd(lambda_n_comp3, v_n_comp3);
                     new_mean_comp3 = mvnrnd(mu_n_comp3',new_covariance_comp3'/k_n_comp3)';
                 end
             end
         end
        
             mean_record_comp1{k,l}{sweep} = new_mean_comp1;
             covariance_record_comp1{k,l}{sweep} = new_covariance_comp1;
             inv_covariance_record_comp1{k,l}{sweep} = inv(new_covariance_comp1);
             
             if(INDICATOR_COMP2)
                 mean_record_comp2{k,l}{sweep} = new_mean_comp2;
                 covariance_record_comp2{k,l}{sweep} = new_covariance_comp2;
                 inv_covariance_record_comp2{k,l}{sweep} = inv(new_covariance_comp2);
             end
             
             if(INDICATOR_COMP3)
                 mean_record_comp3{k,l}{sweep} = new_mean_comp3;
                 covariance_record_comp3{k,l}{sweep} = new_covariance_comp3;
                 inv_covariance_record_comp3{k,l}{sweep} = inv(new_covariance_comp3);
             end
     end;
 end
 
 toc
 
 
  %%%%%%%%%%%%% step 6, sample \alpha and \beta %%%%%%%%%%%%%
  if(ALPHA_BETA_PRIOR==1)
      alpha = gamrnd(a_alpha+K-1,b_alpha-sum(log(1-u(1:K-1))));
      beta = gamrnd(a_beta+K*(L-1),b_beta - sum(sum(log(1-v(:,1:L-1)))));
  end
  
  % record the current parameters values
  zeta(sweep,:)
 
  %%%%%%%%%%%%% step 7, check the figure plot %%%%%%%%%%%%%
  
  if(GRAPHICS && mod(sweep,50)==0)
      
      %plot all subjects - shape component
      figure(11);clf; hold on;
      %plot the shape part;
      color = {'ko', 'ro', 'go', 'bo', 'mo', 'co'};
      hold on;
      labelid=unique(zeta(sweep,:));
      newlable = zeros(NSUB,1);
      for ilab = 1:length(labelid)
          newlable(find(zeta(sweep,:)==labelid(ilab))) = ilab;
      end
      
      for i=1:NSUB
          hold on;
          scatter3(fpca_coeff{i}(1,:),fpca_coeff{i}(2,:),fpca_coeff{i}(3,:),color{min(6,newlable(i))});
      end
      
      %plot all subjects - shape component
      figure(111);clf; hold on;
      for i=1:NSUB
          hold on;
          scatter3(trans{i}(1,:),trans{i}(2,:),trans{i}(3,:),color{min(6,newlable(i))});
      end
      
      
      %plot isub subject - shape component
      figure(12);clf; hold on;
      %plot the shape part;
      isub_id = 1;
      color = {'ko', 'ro', 'go', 'bo', 'mo', 'co','b.'};
      hold on;
      labelid=unique(xi{sweep,isub_id});
      newlabels = xi{sweep,isub_id};
      newlable_isub = zeros(length(newlabels),1);
      for ilab = 1:length(labelid)
          newlable_isub(find(newlabels==labelid(ilab))) = ilab;
      end
      
      for i=1:length(newlabels)
          hold on;
          scatter3(fpca_coeff{isub_id}(1,i),fpca_coeff{isub_id}(2,i),fpca_coeff{isub_id}(3,i),color{min(7,newlable_isub(i))});
      end
      
      figure(13);clf; hold on;
      final_reconstl = input_data{isub_id}.Reconfiber;
      %plot the shape part;
      hold on;
      labelid=unique(xi{sweep,isub_id});
      newlabels = xi{sweep,isub_id};
      newlable_isub = zeros(length(newlabels),1);
      for ilab = 1:length(labelid)
          newlable_isub(find(newlabels==labelid(ilab))) = ilab;
      end
      
      color = {'k', 'r', 'g', 'b', 'm', 'c','k'};
      for i=1:length(newlabels)
          hold on;
          plot3(final_reconstl(1,:,i),final_reconstl(2,:,i),final_reconstl(3,:,i),color{min(7,newlable_isub(i))});
      end
      
      
      %plot another isub subject - shape component
      figure(14);clf; hold on;
      %plot the shape part;
      isub_id = 1;
      color = {'ko', 'ro', 'go', 'bo', 'mo', 'co','b.'};
      hold on;
      labelid=unique(xi{sweep,isub_id});
      newlabels = xi{sweep,isub_id};
      newlable_isub = zeros(length(newlabels),1);
      for ilab = 1:length(labelid)
          newlable_isub(find(newlabels==labelid(ilab))) = ilab;
      end
      
      for i=1:length(newlabels)
          hold on;
          scatter3(fpca_coeff{isub_id}(1,i),fpca_coeff{isub_id}(2,i),fpca_coeff{isub_id}(3,i),color{min(7,newlable_isub(i))});
      end
       
      
  end
end

keyboard;

% MCMC diagnostic
MCMCdiag = 0;
if (MCMCdiag == 1)
    % plot pictures;
    figure, plot(1:length(record),record);
    title('traceplot of k','FontSize',16);
    set(gca,'fontsize',18);
    
    figure, plot(1:length(track_theta),track_theta);
    title('traceplot of \theta','FontSize',16);
    set(gca,'fontsize',18);
    
    % histogram of k
    figure; bar(2:6,ddd(2:6)/sum(ddd));
    set(gca,'fontsize',20);
    title('Histogram of k','FontSize',16);
    
    %autocorrelation plot for \theta
    A_theta = zeros(1,1,TITER);
    A_theta(1,1,:) = track_theta;
    S = mcmcsumm(A_theta);
    mcmctrace(S.acf);
    
    % Effective sample size based on autocorrelation
    disp('Effective sample size of theta:')
    Nh = mcmc_ess_acorr(track_theta)
    
    %Gelman-Rubin R statistic for convergence
    R = mcmcgr(A_theta,10);
    disp('Gelman-Rubin R statistic for convergence:')
    disp(R);
end;