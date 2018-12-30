function [label,means,dist_inClusters] = kmean_curve(dataq,K,initial)

%%%%%%function for K-mean cluster of curves%%%%%%%%%%%
%%%%%%input:  dataq: q function of the data
%%%%%%          K  : cluster numbers
%%%%%       initial: input initilization of the K means

%%%%%output:  means: K-means center
%%%%%         label: labels for each points
%%%%%         dist_inClusters: total distance within each cluster

[dim,M,N] = size(dataq);
t = (0:M-1)/(M-1);
 
q = dataq;
% initial set of k median;
 if(nargin < 3)
    %we need to initialize
    % this is a bad initial
    if(K>1)
        for i=1:K
        initial(:,:,i) = dataq(:,:,ceil(N*rand()));
        end;
    else
        
        % only one clusters;
        [Rend,GAM,means] = mean_curves(q);
        label = ones(1,N);
        dist_inClusters = 0;
        for i= 1:N
            dist_inClusters = dist_inClusters +  pwd_sub(means,q(:,:,i));
        end;
        return;
    end;
 else
     mean = initial;
 end;


mean = initial;

ind = 1;
delta = 100;
ep = 0.1;
tmp_ind = 1;
tmp_dist = zeros(N,K);
while(delta >ep )
    tmp_ind = tmp_ind +1
    if(tmp_ind>14) % Max iteration equal to 15;
        display(delta);
        break;
    end;
    % assign label
    for i=1:N
        for k=1:K
           tmp_dist(i,k) = pwd_sub(mean(:,:,k),q(:,:,i));  
        end;
        ttmp = squeeze(tmp_dist(i,:));
        [dist label(i)] = min(ttmp);
    end;
    
    %check the label;
    if(max(label)~=K)
        display('Only K-1 classes are found;');
    end;
    %find the median;
    for k=1:K
        cluster_data_q = q(:,:,find(label==k));
        [Rend,GAM,mean_updated(:,:,k)]=mean_curves(cluster_data_q);
    end;
    ind = ind +1;
    delta = 0;
    
    for k=1:K
    delta = delta + norm(squeeze(mean_updated(:,:,k)) - squeeze(mean(:,:,k)));
    end;
    
    mean = mean_updated;
end;

means = mean;
%calculate the inside distance
dist_inClusters = zeros(1,K);
for k=1:K
    cluster_idx = find(label==k);
    cluster_data_q = q(:,:,cluster_idx);
    for j=1:length(cluster_idx)
        temp_dist = pwd_sub(mean(:,:,k),cluster_data_q(:,:,j));
        dist_inClusters(1,k) = dist_inClusters(k) + temp_dist;
    end;
end;