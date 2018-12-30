function [Rend,GAM,meanq] = mean_curves(dataq)
%function for caculating the mean of given curvers;
%input: dataq: dim*M*N matrix for N q-function of curves, 
%              where dim is the dimension, M is the length
%              of points along a curve and N is the number of curves. 
% output: Gam: \gamma functions to align each curve to the mean;
%         Rend: rotation matrix: 3 X 3 X N for N curves. 
%         meanq: the mean q function for input data.

%parameters for the whole algorithm
MaxItr = 10; % max iteration number for calculate the mean;
sum_dist = inf; % summry of distance between mean and data;

[dim,M,N] = size(dataq);
t = (0:M-1)/(M-1);
 
%
if N == 0
    display('No inpurt data for mean');
    Rend = 0;
    GAM = 0;
    meanq = 0;
    return;
end;

q=dataq;

% f0 and q0 and original data
q0 = q;

% set initial and update the initial 
mnq = mean(q,3);
for i = 1:N
    dqq(i) = sqrt(sum(sum((q(:,:,i) - mnq).^2)));
end
[tmp, min_ind] = min(dqq);

mq = q(:,:,min_ind);

% update the initial 
parfor k = 1:N
    [q0n(:,:,k),R(:,:,k)] = Find_Best_Rotation(mq, q0(:,:,k)); 
end
Rmean = Mean_Rotation(R);
mq = Rmean'*mq; 

parfor k = 1:N
    norm1 = norm(q0n(:,:,k)); norm2 = norm(mq);
    gam0 = DynamicProgrammingQ(q0n(:,:,k)/norm1,mq/norm2,0,0);
    gam(k,:) = (gam0-gam0(1))/(gam0(end)-gam0(1));  % slight change on scale
end
gamI = SqrtMeanInverse(gam);
gamI_dev = gradient(gamI, 1/(M-1));
mq = interp1(t', mq', (t(end)-t(1)).*gamI' + t(1), 'linear', 'extrap')'.*(ones(dim,1)*sqrt(gamI_dev));


% find the mean
%progressbar
idx = 1;
for r = 1:MaxItr 
%    progressbar(r/MaxItr);
    %%%% Matching Step %%%%
    clear gam gam_dev;
    % use DP to find the optimal warping for each function w.r.t. the mean
    for k = 1:N
        q2new = q0(:,:,k);
        R(:,:,k) = eye(dim);
        for i = 1:5
            [q2n,Rt(:,:,i)] = Find_Best_Rotation(mq(:,:,r),q2new);
            R(:,:,k) = Rt(:,:,i)*R(:,:,k);
            norm1 = norm(q2n); norm2 = norm(mq(:,:,r));
            tic
            gam0 = DynamicProgrammingQ(q2n/norm1,mq(:,:,r)/norm2,0,0);
            totaltime(idx) = toc;
            idx = idx + 1;
            
            gam(k,:) = (gam0-gam0(1))/(gam0(end)-gam0(1));    % normalization
            gam_dev(k,:) = gradient(gam(k,:), 1/(M-1));
            q2new = interp1(t', q2n', (t(end)-t(1)).*gam(k,:)' + t(1), 'linear', 'extrap')'.*(ones(dim,1)*sqrt(gam_dev(k,:)));
            if norm(Rt(:,:,i)-eye(dim)) < 1e-2
                break;
            end
        end
        q(:,:,k) = q2new;
    end
    
    for k = 1:dim
        for n = 1:N
            dd(k,n) = trapz(t, (mq(k,:,r)-q(k,:,n)).^2);
        end
    end
    sum_dist(r+1) = sum(sum(dd)); 

    %%%% Minimization Step %%%%
    % compute the mean of the matched function
    mq(:,:,r+1) = mean(q,3);
    
    qun(r) = norm(mq(:,:,r+1)-mq(:,:,r))/norm(mq(:,:,r));
    %if (ds(r) < ds(r+1) || qun(r) < 1e-2)
    if (qun(r) < 1e-2)
        break;
    end
end



% center the mean
clear q2n R gam;

for k = 1:N
    q2new = squeeze(q0(:,:,k));
    R(:,:,k) = eye(dim);
    GAM(k,:) = linspace(0,1,M);
    for i = 1:5
        [q2n,Rt(:,:,i)] = Find_Best_Rotation(mq(:,:,end),q2new);
        R(:,:,k) = Rt(:,:,i)*R(:,:,k);   % rotation
        
        norm1 = norm(q2n); norm2 = norm(mq(:,:,end));
        gam0 = DynamicProgrammingQ(q2n/norm1,mq(:,:,end)/norm2,0,0); % wrapping
        gamt(i,:) = (gam0-gam0(1))/(gam0(end)-gam0(1));    % normalization
        
        gamt_dev(i,:) = gradient(gamt(i,:), 1/(M-1));
        GAM(k,:) = interp1(t, GAM(k,:), (t(end)-t(1)).*gamt(i,:) + t(1));
        
        q2new = interp1(t', q2n', (t(end)-t(1)).*gamt(i,:)' + t(1), 'linear', 'extrap')'.*(ones(dim,1)*sqrt(gamt_dev(i,:)));
        if norm(Rt(:,:,i)-eye(dim)) < 1e-2
            break;
        end
    end
end

% find the mean element
Rmean = Mean_Rotation(R);    
mqn = Rmean'*mq(:,:,end); % first rotate the mean
GAMI = SqrtMeanInverse(GAM);
GAMI_dev = gradient(GAMI, 1/(M-1));
% second wrap the mean;
mq(:,:,end+1) = interp1(t', mqn', (t(end)-t(1)).*GAMI' + t(1), 'linear', 'extrap')'.*(ones(dim,1)*sqrt(GAMI_dev));


% find the best rotation for 
for k = 1:N
        GAM(k,:) = interp1(t, GAM(k,:), (t(end)-t(1)).*GAMI + t(1));
        R(:,:,k) = Rmean'*R(:,:,k);
end 


%output
meanq = mq(:,:,end);
Rend = R;
