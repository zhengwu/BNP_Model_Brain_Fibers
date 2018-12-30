function [label,aligned_f,dist_inClusters] = kmean_curve_aligned(data,K)

%%%%%%function for K-mean cluster of curves%%%%%%%%%%%
%%%%%%input:  data: original data
%%%%%%          K  : cluster numbers

%%%%%output:  aligned_f: aligned fibers
%%%%%         label: labels for each points
%%%%%         dist_inClusters: total distance within each cluster



% process the data
[dim,M,N] = size(data);
t = linspace(0,1,M);
p = 6;    % number of values to pad 
for i = 1:N;
    f_true(:,:,i) = data(:,:,i);
    f_ori(:,:,i) = [f_true(:,1,i)*ones(1,p) f_true(:,:,i) f_true(:,end,i)*ones(1,p)];
end
t = linspace(0,1,size(f_ori,2));

[dim, M, N] = size(f_ori);

% smooth the data
gk = normpdf(-3:1:3);
gk = gk/sum(gk);
for i = 1:N
    sf(:,:,i) = conv2(f_ori(:,:,i), gk, 'same');
end
clear f;
for i = 1:N
    f(:,:,i) = sf(:,(p+1):end-p,i);
end

[dim, M, N] = size(f);
t = linspace(0,1,M);

% get the SRVF q-function
for i = 1:N
    for k = 1:dim
        q(k,:,i) = gradient(f(k,:,i), 1/(M-1));
    end
    q(:,:,i) = q(:,:,i)./(ones(dim,1)*(sum(q(:,:,i).^2,1).^(1/4))+eps);
end
q_processed = double(q); %SRVF of the processed functions;
f_processed = f; % processed functions;


[label,means_q,dist_inClusters] = kmean_curve(q_processed,K);

% align original functions to the mean function
for k = 1:K
    idx = find(label ==k);
    q_fibers = q_processed(:,:,idx);
    [R,GAM,meanq] = mean_curves(q_fibers);
    
    N_fibers = length(idx);
    
    % align each curve to the mean
    for i = 1:N_fibers
        fn(:,:,idx(i)) = R(:,:,i)*data(:,:,idx(i)); %rotation;
        sl = fn(:,:,idx(i));
        [arclen,~] = arclength(sl(1,:),sl(2,:),sl(3,:));
        %fn(:,:,idx(i)) = fn(:,:,idx(i))/arclen; %scaling
        
        fn(:,:,idx(i)) = interp1(t', fn(:,:,idx(i))', (t(end)-t(1)).*GAM(i,:)' + t(1), 'linear', 'extrap')';
    end
    clear idx fiber;
end;
aligned_f = fn;

%keyboard;

%test Huang's program
%meanf = mean(aligned_f,3);


