function [class_id, K_record, components, alpha_record] = sampler_crp(input_data, num_sweeps)

% Generates samples from the infinite Gaussian mixture model posterior using Chinese Resturant Process

% Input variables are the training_data and the number 

%recover the paramters for each components
INDICATOR_COMP1 = 1; %shape component
INDICATOR_COMP2 = 0; %translation component
INDICATOR_COMP3 = 0; %rotation component
                 
GRAPHICS = 0; % show trace plots and current clustering

fpca_coeff = input_data.fpca_coeff;
trans = input_data.Trans;
rots = input_data.Rots;
mappedrots = input_data.mapedR;
final_reconstl = input_data.Reconfiber;

%pre-processing
trans = trans - mean(trans,2)*ones(1,size(trans,2));
fpca_coeff = fpca_coeff - mean(fpca_coeff,2)*ones(1,size(fpca_coeff,2));
mappedrots = mappedrots - mean(mappedrots,2)*ones(1,size(mappedrots,2));
N = size(fpca_coeff,2);

%component1 = fpca_coeff;
component1 = fpca_coeff;
component2 = trans;
component3 = mappedrots;


if(nargin < 2)
    num_sweeps = 20000;
end



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


% set alpha gamma prior parameters
a_0 = 1;
b_0 = 2;


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
alpha = 0.5;
K_plus = 1;
class_id = zeros(N,num_sweeps);
K_record = zeros(num_sweeps,1);
alpha_record = zeros(num_sweeps,1);
lP_record = zeros(num_sweeps,1);

% seat the first customer at the first table
class_id(:,1) = 1;
% phi{1} = {gaussian(y(:,1),eye(size(y,1)))};
mean_record_comp1{1} = component1(:,1);
covariance_record_comp1{1} = zeros(D1,D1,1);
covariance_record_comp1{1}(:,:,1) = eye(size(component1,1));
inv_covariance_record_comp1{1} = zeros(D1,D1,1);
inv_covariance_record_comp1{1}(:,:,1) = eye(size(component1,1));

mean_record_comp2{1} = component2(:,1);
covariance_record_comp2{1} = zeros(D2,D2,1);
covariance_record_comp2{1}(:,:,1) = eye(size(component2,1));
inv_covariance_record_comp2{1} = zeros(D2,D2,1);
inv_covariance_record_comp2{1}(:,:,1) = eye(size(component2,1));

mean_record_comp3{1} = component3(:,1);
covariance_record_comp3{1} = zeros(D3,D3,1);
covariance_record_comp3{1}(:,:,1) = eye(size(component3,1));
inv_covariance_record_comp3{1} = zeros(D3,D3,1);
inv_covariance_record_comp3{1}(:,:,1) = eye(size(component3,1));


% compute the log likelihood
lp = lp_crp(class_id(:,1),alpha);
for(k=1:K_plus)
    class_k_datapoint_indexes = find(class_id(:,1)==k);
    %     lp= lp + sum(fvnlp(y(:,class_k_datapoint_indexes),...
    %         phi{1}{k}.mean,phi{1}{k}.covariance));
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
    
    class_id(:,sweep) = class_id(:,sweep-1);

    % for each datapoint, unseat it and reseat it, potentially generating a
    % new table
    for(i=1:N)
        m_k = zeros(K_plus,1);

        % compute the CRP prior
        for(k=1:K_plus)
            if(k>K_plus)
                break;
            end
            m_k(k) = length(find(class_id(:,sweep)==k));
            if(class_id(i,sweep)==k)
                m_k(k) = m_k(k)-1;
                if(m_k(k)==0)
                    % delete
                    %                     temp = phi{sweep};
                    %                     phi{sweep} = temp([1:k-1 k+1:end]);

                    mean_record_comp1{sweep} = mean_record_comp1{sweep}(:,[1:k-1 k+1:end]);
                    covariance_record_comp1{sweep} = covariance_record_comp1{sweep}(:,:,[1:k-1 k+1:end]);
                    inv_covariance_record_comp1{sweep} = inv_covariance_record_comp1{sweep}(:,:,[1:k-1 k+1:end]);
                    
                    if(INDICATOR_COMP2)
                        mean_record_comp2{sweep} = mean_record_comp2{sweep}(:,[1:k-1 k+1:end]);
                        covariance_record_comp2{sweep} = covariance_record_comp2{sweep}(:,:,[1:k-1 k+1:end]);
                        inv_covariance_record_comp2{sweep} = inv_covariance_record_comp2{sweep}(:,:,[1:k-1 k+1:end]);
                    end
                    
                    if(INDICATOR_COMP3)
                        mean_record_comp3{sweep} = mean_record_comp3{sweep}(:,[1:k-1 k+1:end]);
                        covariance_record_comp3{sweep} = covariance_record_comp3{sweep}(:,:,[1:k-1 k+1:end]);
                        inv_covariance_record_comp3{sweep} = inv_covariance_record_comp3{sweep}(:,:,[1:k-1 k+1:end]);
                    end
                    
%                    mean_record_comp3{sweep} = mean_record_comp3{sweep}(:,[1:k-1 k+1:end]);
%                    kappa_comp3{sweep} = kappa_comp3{sweep}(:,[1:k-1 k+1:end]);
                    
                    for(j=k:K_plus)
                        change_inds = find(class_id(:,sweep)==j);
                        class_id(change_inds,sweep) = j-1;
                    end
                    K_plus = K_plus-1;
                    m_k = m_k(1:end-1);
                end
            end
        end

        % sneakily add on the new table generation prior prob to the vec.
        prior = [m_k; alpha]/(N-1+alpha);
        likelihood = zeros(length(prior)-1,1);
        

        % compute the per class data likelihood p(x_i|g_i,\Theta)
        for(l = 1:length(likelihood))

            %             likelihood(l) = lp(temp{l},y(:,i));
            %                 [p lp] = p(temp{l},y(:,i));
            %                 lp = fvnlp(y(:,1),temp{l}.mean, temp{l}.covariance);
            %                 lp = fvnlp(y(:,i),mean_record{sweep}(:,l),covariance_record_comp1{sweep}(:,:,l));
            l1 = component1(:,i)-mean_record_comp1{sweep}(:,l);
            lp1 = -l1'*inv_covariance_record_comp1{sweep}(:,:,l)*l1/2-.5*log(det(covariance_record_comp1{sweep}(:,:,l)))-(D1/2)*log(2*pi);

            if(INDICATOR_COMP2)
                l2 = component2(:,i)-mean_record_comp2{sweep}(:,l);
                lp2 = -l2'*inv_covariance_record_comp2{sweep}(:,:,l)*l2/2-.5*log(det(covariance_record_comp2{sweep}(:,:,l)))-(D2/2)*log(2*pi);
            end
            
            if(INDICATOR_COMP3)
                l3 = component3(:,i) - mean_record_comp3{sweep}(:,l);
                lp3 = -l3'*inv_covariance_record_comp3{sweep}(:,:,l)*l3/2-.5*log(det(covariance_record_comp3{sweep}(:,:,l)))-(D3/2)*log(2*pi);
            end
                       
            if(INDICATOR_COMP3&&INDICATOR_COMP2)
                likelihood(l) = lp1+lp2+lp3;
            elseif(INDICATOR_COMP2)
                likelihood(l) = lp1+lp2;
            elseif(INDICATOR_COMP3)
                likelihood(l) = lp1+lp3;
            else 
                likelihood(l) = lp1;
            end
            
        end

     
        new_table_likelihood = 0;
        
        new_table_covariance1 = iwishrnd(lambda_0_comp1, v_0_comp1);
        new_table_mean1 = mvnrnd(mu_0_comp1',new_table_covariance1'/k_0_comp1)';
        l1 = component1(:,i)-new_table_mean1;
        inverse_new_table_covariance1 = inv(new_table_covariance1);
        new_table_likelihood = new_table_likelihood-l1'*inverse_new_table_covariance1*l1/2-.5*log(det(new_table_covariance1))-(D1/2)*log(2*pi);
        
        % for component2
        if(INDICATOR_COMP2)
            new_table_covariance2 = iwishrnd(lambda_0_comp2, v_0_comp2);
            new_table_mean2 = mvnrnd(mu_0_comp2',new_table_covariance2'/k_0_comp2)';
            l2 = component2(:,i)-new_table_mean2;
            inverse_new_table_covariance2 = inv(new_table_covariance2);
            new_table_likelihood = new_table_likelihood-l2'*inverse_new_table_covariance2*l2/2-.5*log(det(new_table_covariance2))-(D2/2)*log(2*pi);
        end
        
        %for component3 
        if(INDICATOR_COMP3)
            new_table_covariance3 = iwishrnd(lambda_0_comp3, v_0_comp3);
            new_table_mean3 = mvnrnd(mu_0_comp3',new_table_covariance3'/k_0_comp3)';
            l3 = component3(:,i)-new_table_mean3;
            inverse_new_table_covariance3 = inv(new_table_covariance3);
            new_table_likelihood = new_table_likelihood-l3'*inverse_new_table_covariance3*l3/2-.5*log(det(new_table_covariance3))-(D3/2)*log(2*pi);
        end
         
        likelihood = [likelihood; new_table_likelihood];
        likelihood = exp(likelihood-max(likelihood));%/sum(exp(likelihood-max(likelihood)));
        likelihood = likelihood/sum(likelihood);
        if(sum(likelihood==0))
            likelihood(end) = eps;
        end
        
        % compute the posterior over seating assignment for data i
        temp = prior.*likelihood;
        temp = temp/sum(temp);

        cdf = cumsum(temp);
        rn = rand;
        new_class_id = min(find((cdf>rn)==1));
        %         picker = multinomial(temp);
        %         new_class_id = sample(picker,1);

        % if the new table was chosen, add it to the list of tables a
        % datum can sit at
        %         temp = phi{sweep};
        if(new_class_id == K_plus+1)
            %             phi{sweep} = {temp{:}, new_table_density};
            mean_record_comp1{sweep} = [mean_record_comp1{sweep} new_table_mean1];
            covariance_record_comp1{sweep}(:,:,new_class_id) =  new_table_covariance1;
            inv_covariance_record_comp1{sweep}(:,:,new_class_id) =  inv(new_table_covariance1);
            
            if(INDICATOR_COMP2)
                mean_record_comp2{sweep} = [mean_record_comp2{sweep} new_table_mean2];
                covariance_record_comp2{sweep}(:,:,new_class_id) =  new_table_covariance2;
                inv_covariance_record_comp2{sweep}(:,:,new_class_id) =  inv(new_table_covariance2);
            end
            
            if(INDICATOR_COMP3)
                mean_record_comp3{sweep} = [mean_record_comp3{sweep} new_table_mean3];
                covariance_record_comp3{sweep}(:,:,new_class_id) =  new_table_covariance3;
                inv_covariance_record_comp3{sweep}(:,:,new_class_id) =  inv(new_table_covariance3);
            end
            
            K_plus = K_plus+1;
        end

        % record the new table
        class_id(i,sweep) = new_class_id;

    end

    % flag to set whether to use a normal inverse wishart prior
    % or a Jeffries non-informative prior...  the Jeffries prior is not
    % tested
    NIW = 1;

    % update the cluster parameters
    for(k=1:K_plus)
        class_k_datapoint_indexes = find(class_id(:,sweep)==k);
        N_this_class = length(class_k_datapoint_indexes);
        
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


    SAMPLEALPHA = 0;
    METROPOLISALPHA = 1;
    
    if(~SAMPLEALPHA)
        % compute the log likelihood
        lp = lp_crp(class_id(:,sweep),alpha);
        for(k=1:K_plus)
            class_k_datapoint_indexes = find(class_id(:,sweep)==k);
            %         lp= lp + sum(fvnlp(y(:,class_k_datapoint_indexes),phi{sweep}{k}.mean,phi{sweep}{k}.covariance));
            lp= lp + sum(fvnlp(component1(:,class_k_datapoint_indexes),mean_record_comp1{sweep}(:,k),covariance_record_comp1{sweep}(:,:,k)));
            
            if(INDICATOR_COMP2)
                lp= lp + sum(fvnlp(component2(:,class_k_datapoint_indexes),mean_record_comp2{sweep}(:,k),covariance_record_comp2{sweep}(:,:,k)));
            end
            
            if(INDICATOR_COMP3)
                lp= lp + sum(fvnlp(component3(:,class_k_datapoint_indexes),mean_record_comp3{sweep}(:,k),covariance_record_comp3{sweep}(:,:,k)));
            end
            %         mu = phi{sweep}{k}.mean;
            %         sigma = phi{sweep}{k}.covariance;
            mu = mean_record_comp1{sweep}(:,k);
            sigma = covariance_record_comp1{sweep}(:,:,k);
            lp = lp + lpnormalinvwish(mu,sigma,mu_0_comp1,k_0_comp1,v_0_comp1,lambda_0_comp1);
            
            if(INDICATOR_COMP2)
                mu = mean_record_comp2{sweep}(:,k);
                sigma = covariance_record_comp2{sweep}(:,:,k);
                lp = lp + lpnormalinvwish(mu,sigma,mu_0_comp2,k_0_comp2,v_0_comp2,lambda_0_comp2);
            end
            
            if(INDICATOR_COMP3)
                mu = mean_record_comp3{sweep}(:,k);
                sigma = covariance_record_comp3{sweep}(:,:,k);
                lp = lp + lpnormalinvwish(mu,sigma,mu_0_comp3,k_0_comp3,v_0_comp3,lambda_0_comp3);
            end    
        end


    end
    % I could not get the Gibbs sampler for alpha to work.
    if(~METROPOLISALPHA && SAMPLEALPHA == 1)
    %sample alpha
        nu = betarnd(alpha+1,N);
        a = 1;
        b = 1;
        % this is the same as eqn. 14 of Escobar and West 1994 Bayesian
        % Density Estimation and Inference Using Mixtures
        pia = (a+K_plus-1)/((a+K_plus-1)+(b-log(nu))*N)
        if(rand < pia)
            alpha = gamrnd(a+K_plus,b-log(nu))
        else
            alpha = gamrnd(a+K_plus-1,b-log(nu))
        end
    end

    if(METROPOLISALPHA && SAMPLEALPHA==1)
        % use a metropolis step to sample alpha
        lp = lp_crp(class_id(:,sweep),alpha);
        for(k=1:K_plus)
            class_k_datapoint_indexes = find(class_id(:,sweep)==k);
            %         lp= lp + sum(fvnlp(y(:,class_k_datapoint_indexes),phi{sweep}{k}.mean,phi{sweep}{k}.covariance));
            lp= lp + sum(fvnlp(component1(:,class_k_datapoint_indexes),mean_record_comp1{sweep}(:,k),covariance_record_comp1{sweep}(:,:,k)));
            
            if(INDICATOR_COMP2)
                lp= lp + sum(fvnlp(component2(:,class_k_datapoint_indexes),mean_record_comp2{sweep}(:,k),covariance_record_comp2{sweep}(:,:,k)));
            end
            
            if(INDICATOR_COMP3)
                lp= lp + sum(fvnlp(component3(:,class_k_datapoint_indexes),mean_record_comp3{sweep}(:,k),covariance_record_comp3{sweep}(:,:,k)));
            end
            %         mu = phi{sweep}{k}.mean;
            %         sigma = phi{sweep}{k}.covariance;
            mu = mean_record_comp1{sweep}(:,k);
            sigma = covariance_record_comp1{sweep}(:,:,k);
            lp = lp + lpnormalinvwish(mu,sigma,mu_0_comp1,k_0_comp1,v_0_comp1,lambda_0_comp1);
            
            if(INDICATOR_COMP2)
            mu = mean_record_comp2{sweep}(:,k);
            sigma = covariance_record_comp2{sweep}(:,:,k);
            lp = lp + lpnormalinvwish(mu,sigma,mu_0_comp2,k_0_comp2,v_0_comp2,lambda_0_comp2);
            end
            
            if(INDICATOR_COMP3)
                mu = mean_record_comp3{sweep}(:,k);
                sigma = covariance_record_comp3{sweep}(:,:,k);
                lp = lp + lpnormalinvwish(mu,sigma,mu_0_comp3,k_0_comp3,v_0_comp3,lambda_0_comp3);
            end
            
            lp = lp + -  gamlike([a_0 b_0],alpha);
        end

        alpha_prop = alpha + randn*.1;

        if(alpha_prop < 0)
            lp_alpha_prop = -Inf;
        else
            lp_alpha_prop = lp_crp(class_id(:,sweep),alpha_prop);
            lp_alpha_prop = lp_alpha_prop -  gamlike([a_0 b_0],alpha_prop);
        end
        log_acceptance_ratio = lp_alpha_prop - lp;
        if(log(rand)<min(log_acceptance_ratio,0))
            alpha = alpha_prop;
        end
    end

    % record the current parameters values
    K_record(sweep) = K_plus;
    alpha_record(sweep) = alpha;
    lP_record(sweep) = lp;
    if(GRAPHICS&mod(sweep,10)==0) % every 100 iterations, display the figure
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


        figure(1)
        %         plot_mixture(component1(1:2,1:32:end),class_id(1:32:end,sweep))
        plot_mixture3(component1(1:3,:),class_id(:,sweep))
        
        if(INDICATOR_COMP2)
        figure(2)
        %         plot_mixture(component1(1:2,1:32:end),class_id(1:32:end,sweep))
        plot_mixture3(component2(1:3,:),class_id(:,sweep))
        end
        drawnow
        
        if(INDICATOR_COMP3)
            figure(3)
            plot_mixture3(component3(1:3,:),class_id(:,sweep))
        end
        %     pause(1)
        
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
