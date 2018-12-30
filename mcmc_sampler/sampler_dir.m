function [class_id, K_record, components, alpha_record] = sampler_dir(input_data, K, num_sweeps)

% Generates samples from the Gaussian mixture model posterior

% Input variables are training_data, initial K, and the number MCMC samples
% Output: class_id: posterior samples of memberships
%         K_record: posterior samples of # of clusters
%         components: posterior samples of mean and covariance of each components 
%         alpha_record: hyperparameter alpha

fpca_coeff = input_data.fpca_coeff;
trans = input_data.Trans;
mappedrots = input_data.mapedR;
final_reconstl = input_data.Reconfiber;

trans = trans - mean(trans,2)*ones(1,size(trans,2));
fpca_coeff = fpca_coeff - mean(fpca_coeff,2)*ones(1,size(fpca_coeff,2));
mappedrots = mappedrots - mean(mappedrots,2)*ones(1,size(mappedrots,2));
N = size(fpca_coeff,2);

%recover the paramters for each components
INDICATOR_COMP1 = 1; 
INDICATOR_COMP2 = 1;
INDICATOR_COMP3 = 1;
                
component1 = fpca_coeff;
component2 = trans;
component3 = mappedrots;


if(nargin < 3)
    num_sweeps = 20000;
elseif(nargin < 2)
    K = 10;
    num_sweeps = 20000;
end

GRAPHICS = 0; % show trace plots and current clustering


%%%%%%%%%%%%%%%% set the hyper parameters%%%%%%%%%%%%%%%%
% set normal inverse wishart hyper parameters and gamma prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set normal inverse wishart hyper parameters for comp1
mu_0_comp1 = zeros(size(component1(:,1)));
k_0_comp1 = 1;
v_0_comp1 = size(component1(:,1),1);
lambda_0_comp1 = eye(size(component1(:,1),1))*0.5;


% set normal inverse wishart hyper parameters for comp2
mu_0_comp2 = zeros(size(component2(:,1)));
k_0_comp2 = 1;
v_0_comp2 = size(component2(:,1),1);
lambda_0_comp2 = eye(size(component2(:,1),1))*0.5;


% set hyper parameters for comp3
mu_0_comp3 = zeros(size(component3(:,1)));
k_0_comp3 = 1;
v_0_comp3 = size(component3(:,1),1);
lambda_0_comp3 = eye(size(component3(:,1),1))*0.5;



%%%%%%%%%%%%%%%% initialization %%%%%%%%%%%%%%%%
% initalize the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for component1
D1 = size(component1,1);
mean_record_comp1 = cell(num_sweeps,1);
covariance_record_comp1 = cell(num_sweeps,1);
inv_covariance_record_comp1 = cell(num_sweeps,1);
lP_record_comp1 = zeros(num_sweeps,1);

%for component2
D2 = size(component2,1);
mean_record_comp2 = cell(num_sweeps,1);
covariance_record_comp2 = cell(num_sweeps,1);
inv_covariance_record_comp2 = cell(num_sweeps,1);
lP_record_comp2 = zeros(num_sweeps,1);

%for component3
D3 = size(component3,1);
mean_record_comp3 = cell(num_sweeps,1);
covariance_record_comp3 = cell(num_sweeps,1);
inv_covariance_record_comp3 = cell(num_sweeps,1);
lP_record_comp3 = zeros(num_sweeps,1);

%shared in all components
N = size(component1,2);
pi_sample = zeros(num_sweeps,K);
for i = 1:K
    pi_sample(1,i) = 1/K;
end
alpha = 0.01*K;

class_id = zeros(N,num_sweeps);
noempty_K_record = zeros(num_sweeps,1);
alpha_record = zeros(num_sweeps,1);
lP_record = zeros(num_sweeps,1);

%%%%%%%%%%%%%%%% initialize parameters%%%%%%%%%%%%%%%%
% initialize index
for i = 1:N
    class_id(i,1) = 1+mod(i,2);
end
class_id(:,1) =  class_id(randperm(N),1);
 
%for each center, initialize its parameters
for k=1:K
    fiberK_ind = find(class_id(:,1)==k);
    if(length(fiberK_ind)==0)
        %assign priors
            %for component 1
            new_covariance_comp1 = iwishrnd(inv(lambda_0_comp1), v_0_comp1+2);
            mean_record_comp1{1}(:,k) = mvnrnd(mu_0_comp1',new_covariance_comp1'/k_0_comp1)';%component1{1+mod(k,NSUB)}(:,l);
            covariance_record_comp1{1}(:,:,k) = new_covariance_comp1;
            inv_covariance_record_comp1{1}(:,:,k) = inv(covariance_record_comp1{1}(:,:,k));
            
            %for component 2
            new_covariance_comp2 = iwishrnd(inv(lambda_0_comp2), v_0_comp2+2);
            mean_record_comp2{1}(:,k) = mvnrnd(mu_0_comp2',new_covariance_comp2'/k_0_comp2)';%component1{1+mod(k,NSUB)}(:,l);
            covariance_record_comp2{1}(:,:,k) = new_covariance_comp2;
            inv_covariance_record_comp2{1}(:,:,k) = inv(covariance_record_comp2{1}(:,:,k));
            
            %for component 3
            new_covariance_comp3 = iwishrnd(inv(lambda_0_comp3), v_0_comp3+2);
            mean_record_comp3{1}(:,k) = mvnrnd(mu_0_comp3',new_covariance_comp3'/k_0_comp3)';%component1{1+mod(k,NSUB)}(:,l);
            covariance_record_comp3{1}(:,:,k) = new_covariance_comp3;
            inv_covariance_record_comp3{1}(:,:,k) = inv(covariance_record_comp3{1}(:,:,k));
            
            
    else
            
            kl_comp1 = component1(:,fiberK_ind);
            kl_comp2 = component2(:,fiberK_ind);
            kl_comp3 = component3(:,fiberK_ind);
            
            mean_record_comp1{1}(:,k) = mean(kl_comp1,2);
            covariance_record_comp1{1}(:,:,k) = cov(kl_comp1') + 0.1*eye(3,3);
            inv_covariance_record_comp1{1}(:,:,k) = inv(covariance_record_comp1{1}(:,:,k));
            
            mean_record_comp2{1}(:,k) = mean(kl_comp2,2);
            covariance_record_comp2{1}(:,:,k) = cov(kl_comp2') + 0.1*eye(3,3);
            inv_covariance_record_comp2{1}(:,:,k) = inv(covariance_record_comp2{1}(:,:,k));
            
            mean_record_comp3{1}(:,k) = mean(kl_comp3,2);
            covariance_record_comp3{1}(:,:,k) = cov(kl_comp3') + 0.1*eye(3,3);
            inv_covariance_record_comp3{1}(:,:,k) = inv(covariance_record_comp3{1}(:,:,k));
            
    end
end


% compute the log likelihood
lp = lp_dirichlet(pi_sample(1,:)',(ones(K,1)*alpha/K));
%lp = 0;
for(k=1:K)
    class_k_datapoint_indexes = find(class_id(:,1)==k);
    if(length(class_k_datapoint_indexes)>0)
        lp= lp + sum(fvnlp(component1(:,class_k_datapoint_indexes),...
            mean_record_comp1{1}(:,k),covariance_record_comp1{1}(:,:,k)));
        
        if(INDICATOR_COMP2)
            lp= lp + sum(fvnlp(component2(:,class_k_datapoint_indexes),...
                mean_record_comp2{1}(:,k),covariance_record_comp2{1}(:,:,k)));
        end
        
        if(INDICATOR_COMP3)
            %component3 - rotations
            lp= lp + sum(fvnlp(component3(:,class_k_datapoint_indexes),...
                mean_record_comp3{1}(:,k),covariance_record_comp3{1}(:,:,k)));
            %lp = lp+ LogLikelihoodDataTerm(mean_record_comp3{1}(:,k),kappa_comp3{1}(k),component3(:,:,class_k_datapoint_indexes))
        end
        
        mu = mean_record_comp1{1}(:,k);
        sigma = covariance_record_comp1{1}(:,:,k);
        lp = lp + lpnormalinvwish(mu,sigma,mu_0_comp1,k_0_comp1,v_0_comp1,lambda_0_comp1);
        
        if(INDICATOR_COMP2)
            mu = mean_record_comp2{1}(:,k);
            sigma = covariance_record_comp2{1}(:,:,k);
            lp = lp + lpnormalinvwish(mu,sigma,mu_0_comp2,k_0_comp2,v_0_comp2,lambda_0_comp2);
        end
        
        if(INDICATOR_COMP3)
            mu = mean_record_comp3{1}(:,k);
            sigma = covariance_record_comp3{1}(:,:,k);
            lp = lp + lpnormalinvwish(mu,sigma,mu_0_comp3,k_0_comp3,v_0_comp3,lambda_0_comp3);
        end
    end
    
end
lP_record(1) = lp;


% run the Gibbs sampler
for(sweep = 2:num_sweeps)
    if(mod(sweep,100)==0)
        disp(['Sweep ' num2str(sweep) '/' num2str(num_sweeps)])
        K_plus
    end
    
    mean_record_comp1{sweep} = mean_record_comp1{sweep-1};
    covariance_record_comp1{sweep} = covariance_record_comp1{sweep-1};
    inv_covariance_record_comp1{sweep} = inv_covariance_record_comp1{sweep-1};
    
    if(INDICATOR_COMP2)
        mean_record_comp2{sweep} = mean_record_comp2{sweep-1};
        covariance_record_comp2{sweep} = covariance_record_comp2{sweep-1};
        inv_covariance_record_comp2{sweep} = inv_covariance_record_comp2{sweep-1};
    end
    
    if(INDICATOR_COMP3)
        mean_record_comp3{sweep} = mean_record_comp3{sweep-1};
        covariance_record_comp3{sweep} = covariance_record_comp3{sweep-1};
        inv_covariance_record_comp3{sweep} = inv_covariance_record_comp3{sweep-1};
    end
    
    %class_id(:,sweep) = class_id(:,sweep-1);
    pi_sample(sweep,:) = pi_sample(sweep-1,:);
    
    % step 1 - reassign cluster for each fiber
    for ifb = 1:N
        %for each subject, reassign a cluster to it
        
        %get the data;
        curr_datap_c1 = component1(:,ifb);
        curr_datap_c2 = component2(:,ifb);
        curr_datap_c3 = component3(:,ifb);
        
        %calculate the probability assigning each fiber to a class
        data_likelihood_lp = zeros(K,1);
        for k=1:K
            %prior
            tmp_pi = pi_sample(sweep,k);
            
            %data
            l1 = curr_datap_c1 - mean_record_comp1{sweep}(:,k);
            if(INDICATOR_COMP2)
                l2 = curr_datap_c2 - mean_record_comp2{sweep}(:,k);
            end
            if(INDICATOR_COMP3)
                l3 = curr_datap_c3 - mean_record_comp3{sweep}(:,k);
            end
            
            %likelihood for each component;
            lp1 = -l1'*inv_covariance_record_comp1{sweep}(:,:,k)*l1/2-.5*log(det(covariance_record_comp1{sweep}(:,:,k)))-(D1/2)*log(2*pi);
            if(INDICATOR_COMP2)
                lp2 = -l2'*inv_covariance_record_comp2{sweep}(:,:,k)*l2/2-.5*log(det(covariance_record_comp2{sweep}(:,:,k)))-(D2/2)*log(2*pi);
            end
            
            if(INDICATOR_COMP3)
                %lp3 = LogLikelihoodDataTerm(mean_record_comp3{k,l}{sweep},kappa_comp3{sweep}(l),l3);
                lp3 = -l3'*inv_covariance_record_comp3{sweep}(:,:,k)*l3/2-.5*log(det(covariance_record_comp3{sweep}(:,:,k)))-(D3/2)*log(2*pi);
            end
            
            
            if(INDICATOR_COMP3&&INDICATOR_COMP2)
                data_likelihood_kl_lp(k) = lp1+lp2+lp3;
            elseif(INDICATOR_COMP2)
                data_likelihood_kl_lp(k) = lp1+lp2;
            else
                data_likelihood_kl_lp(k) = lp1;
            end
            
        end
        
        data_likelihood_kl_lp = data_likelihood_kl_lp + log(pi_sample(sweep,:));
       
        likelihood = exp(data_likelihood_kl_lp-max(data_likelihood_kl_lp));%/sum(exp(likelihood-max(likelihood)));
        likelihood = likelihood/sum(likelihood);
        if(sum(likelihood)==0)
            likelihood(end) = eps;
            keyboard;
        end
        
        cdf = cumsum(likelihood);
        rn = rand;
        new_class_id = min(find((cdf>rn)==1));
        
        % record the new clustering id;
        class_id(ifb,sweep) = new_class_id;
    end
   
    % step 2 - update pi
    m = zeros(K,1);
    for k=1:K
        class_k_datapoint_indexes = find(class_id(:,sweep)==k);
        N_this_class = length(class_k_datapoint_indexes);
        m(k) = N_this_class;
    end
    tmp_a = alpha*ones(K,1)/K + m;
    new_pi = drchrnd(tmp_a',1);
    pi_sample(sweep,:) = new_pi;
        
    % flag to set whether to use a normal inverse wishart prior
    % or a Jeffries non-informative prior...  the Jeffries prior is not
    % tested
    NIW = 1;
    
    % step 3 - update the cluster parameters
    for(k=1:K)
        class_k_datapoint_indexes = find(class_id(:,sweep)==k);
        N_this_class = length(class_k_datapoint_indexes);
        m(k) = N_this_class;
        
        if(N_this_class==0) % if there is no elements in a class, we sample from the prior distribution
             new_covariance_comp1 = iwishrnd(lambda_0_comp1, v_0_comp1);
             new_mean_comp1 = mvnrnd(mu_0_comp1',new_covariance_comp1'/k_0_comp1)';
             
             new_covariance_comp2 = iwishrnd(lambda_0_comp2, v_0_comp2);
             new_mean_comp2 = mvnrnd(mu_0_comp2',new_covariance_comp2'/k_0_comp2)';

             new_covariance_comp3 = iwishrnd(lambda_0_comp3, v_0_comp3);
             new_mean_comp3 = mvnrnd(mu_0_comp3',new_covariance_comp3'/k_0_comp3)';
        else
             
            S1 = cov(component1(:,class_k_datapoint_indexes)')'*(N_this_class-1);
            y_bar1 = mean(component1(:,class_k_datapoint_indexes)')';
            mu_n_comp1 = (N_this_class*y_bar1+...
                k_0_comp1*mu_0_comp1)/(k_0_comp1+N_this_class);
        
            if(INDICATOR_COMP2)
                S2 = cov(component2(:,class_k_datapoint_indexes)')'*(N_this_class-1);
                y_bar2 = mean(component2(:,class_k_datapoint_indexes)')';
                mu_n_comp2 = (N_this_class*y_bar2+...
                    k_0_comp2*mu_0_comp2)/(k_0_comp2+N_this_class);
            end
            
            if(INDICATOR_COMP3)
                S3 = cov(component3(:,class_k_datapoint_indexes)')'*(N_this_class-1);
                y_bar3 = mean(component3(:,class_k_datapoint_indexes)')';
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
            else
                new_covariance = iwishrnd(S,N_this_class-1);
                new_mean = mvnrnd(mean(component1(:,class_k_datapoint_indexes)')',...
                    new_covariance'/N_this_class)';
            end
            mean_record_comp1{sweep}(:,k) = new_mean_comp1;
            covariance_record_comp1{sweep}(:,:,k) = new_covariance_comp1;
            inv_covariance_record_comp1{sweep}(:,:,k) = inv(new_covariance_comp1);
            
            
            if(INDICATOR_COMP2)
                mean_record_comp2{sweep}(:,k) = new_mean_comp2;
                covariance_record_comp2{sweep}(:,:,k) = new_covariance_comp2;
                inv_covariance_record_comp2{sweep}(:,:,k) = inv(new_covariance_comp2);
            end
            
            if(INDICATOR_COMP3)
                mean_record_comp3{sweep}(:,k) = new_mean_comp3;
                covariance_record_comp3{sweep}(:,:,k) = new_covariance_comp3;
                inv_covariance_record_comp3{sweep}(:,:,k) = inv(new_covariance_comp3);
            end
        end
    end

    
    % record the current parameters values
    K_plus = length(unique(class_id(:,sweep)));
    K_record(sweep) = K_plus;
    alpha_record(sweep) = alpha;
    lP_record(sweep) = lp;
    if(GRAPHICS&mod(sweep,100)==0) % every 100 iterations, display the figure
        figure(4)
        subplot(3,1,1)
        plot(1:sweep,lP_record(1:sweep));
        title('Log P')
        subplot(3,1,2)
        plot(1:sweep,K_record(1:sweep));
        title('K');
        subplot(3,1,3)
        plot(1:sweep,alpha_record(1:sweep));
        title('alpha');


        figure(1);clf;
        plot_mixture3(component1(1:3,:),class_id(:,sweep))
        
        if(INDICATOR_COMP2)
        figure(2);clf;
        plot_mixture3(component2(1:3,:),class_id(:,sweep))
        end
        drawnow
        
        if(INDICATOR_COMP3)
            figure(3);clf;
            plot_mixture3(component3(1:3,:),class_id(:,sweep))
        end
        
        figure(10);clf;
        newlable_isub = class_id(:,sweep);
        color = {'k', 'r', 'g', 'b', 'm', 'c'};
        for i=1:N
            hold on;
            plot3(final_reconstl(1,:,i),final_reconstl(2,:,i),final_reconstl(3,:,i),color{min(6,newlable_isub(i))});
        end
        
        
    end


end

components.postermean_c1 = mean_record_comp1;
components.postervar_c1 = covariance_record_comp1;

components.postermean_c2 = mean_record_comp2;
components.postervar_c2 = covariance_record_comp2;

components.postermean_c3 = mean_record_comp3;
components.postervar_c3 = covariance_record_comp3;



%%%%%%%%%%%%%%%%%%%%%% MCMC diagnosis
MCMCdiag = 0;
if(MCMCdiag == 1)
    %post-processing to fix the lable switching issue;
    record_label = class_id';
    burnin = 1000;
    N = size(record_label,2);
    
    %plot the original data;
    finalB = zeros(N,N);
    for i=burnin:num_sweeps
        tempB = calculate_B(squeeze(record_label(i,:)),N);
        finalB = finalB + tempB;
        psoteriorK(i-burnin+1) = length(unique(record_label(i,:)));
    end
    
    [c,cn,thrd] = BtoCluster_zw(finalB,N,mode(psoteriorK));
    
    %for the first cluster, pull out its data from the chain;
    indx1 = find(c==1);
    all_sample_ind1 = record_label(:,indx1);
    for i=1:size(all_sample_ind1,1)
        ind_cluster1(i) =  mode(all_sample_ind1(i,:));
        pi_cluster1(i) = pi_sample(i,ind_cluster1(i));
        mean_cluster1(:,i) = mean_record_comp1{i}(:,ind_cluster1(i));
    end
    
    %for the second cluster, pull out its data from the chain;
    indx2 = find(c==2);
    all_sample_ind2 = record_label(:,indx2);
    for i=1:size(all_sample_ind2,1)
        ind_cluster2(i) =  mode(all_sample_ind2(i,:));
        pi_cluster2(i) = pi_sample(i,ind_cluster2(i));
        mean_cluster2(:,i) = mean_record_comp1{i}(:,ind_cluster2(i));
    end
    
    %MCMC diagnosis
    for i=1000:num_sweeps
        mean_record_comp1_x(i) = mean_record_comp1{i}(1,1);
    end
    
    % save result for Gelman-Rubin R
    %save run5 ind_cluster1 pi_cluster1 mean_cluster1 ind_cluster2 pi_cluster2 mean_cluster2;
    
    
    % plot pictures;
    figure, plot(1:length(K_record(burnin+1:end)),K_record(burnin+1:end),'linewidth',2);
    title('Traceplot of K','FontSize',16);
    set(gca,'fontsize',18);
    
    figure, plot(1:length(pi_cluster1(burnin+1:end)),pi_cluster1(burnin+1:end));
    hold on, plot(1:length(pi_cluster2(burnin+1:end)),pi_cluster2(burnin+1:end),'r');
    title('Traceplot of Pi','FontSize',16);
    set(gca,'fontsize',18);
    
    
    figure, plot(1:num_sweeps-burnin,mean_cluster1(:,burnin+1:end));
    title('traceplot of mean 1','FontSize',16);
    set(gca,'fontsize',18);
    
    figure, plot(1:num_sweeps-burnin,mean_cluster2(:,burnin+1:end));
    title('traceplot of mean 2','FontSize',16);
    set(gca,'fontsize',18);
    
    %effect sample size;
    disp('Effective sample size of theta:')
    Nh = mcmc_ess_acorr(pi_cluster1(burnin+1:end));
    
    disp('Effective sample size of theta:')
    Nh = mcmc_ess_acorr(pi_cluster2(burnin+1:end));
    
    Nh = mcmc_ess_acorr(mean_cluster1(1,burnin+1:end));
    
    
    %     %Gelman-Rubin R statistic for convergence
    %     R = mcmcgr(A_theta,10);
    %     disp('Gelman-Rubin R statistic for convergence:')
    %     disp(R);

    
end

