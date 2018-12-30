function fdiff = pca_analysis(fibers,all_len)
% function for performing pca analysis
% input: fibers - 3*M*N matrix, representing N fibers, each fiber
% is 3 dimensional with M points;
%        all_len - fiber lengths

% output: fdiff - vector of 30, representing mean difference between N pc
% components recovered fibers and orginal fibers

%display aligned fibers
rcentered_stl = fibers;
[dim,M,N]=size(rcentered_stl);


%plot 3D fibers
figure(10);clf; hold on;
for i=1:N
    plot3(rcentered_stl(1,:,i),rcentered_stl(2,:,i),rcentered_stl(3,:,i),'linewidth',3);
end
xlabel('X axis');ylabel('Y axis');zlabel('Z axis');
set(gca,'fontsize',16);
axis off;
view([164,22])


%plot x-y-z components
cl = {[1 0.8 0], [0.3 0.3 0.6], [1 0.5 0.5], 'm', 'c', [0,0.4,0.5], 'g', 'b', 'r','y', 'k'};
figure(13);clf;
for k = 1:3
    subplot(3,1,k);
    for n = 1:N
        plot(1:size(rcentered_stl(:,:,1),2), rcentered_stl(k,:,n), 'b','linewidth',0.8); hold on;
    end
    
    hold on;
    
    if(k==1)
        title(sprintf('dim-x'), 'fontsize', 18);
    elseif(k==2)
        title(sprintf('dim-y'), 'fontsize', 18);

    elseif(k==3)
        title(sprintf('dim-z'), 'fontsize', 18);
    end
    
    axis tight;
    set(gca, 'fontsize', 24);
end


% %%%%%%%%%%%%%%%%% pca to generate basis %%%%%%%%%%%%%%%%%%%%
%parameters for pca
nocomp = 30; % number of principal components; right now, it is fixed; 
%adaptive nocomp can be used. 

NP = 1:nocomp;

    
%parameters to save pca result;
K = 1;
pca_basis = zeros(dim,M,nocomp); %all subjects share the same basis
pca_engvalues = zeros(dim,M);

%get the data ready
currt_mn = mean(rcentered_stl,3);

%remove mean;
for i=1:N
    normalized_f(:,:,i) = rcentered_stl(:,:,i) - currt_mn;
end
    
%do pca for each component
for cp = 1:3
    cp_f = squeeze(normalized_f(cp,:,:))';
    
    K = cov(cp_f);
    [U,S,~] = svd(K);
    sdiag = diag(S);
    stdS = sqrt(sdiag);
    
    % PCA basis
    pca_basis(cp,:,:) = U(:,NP);
    pca_engvalues(cp,:) = sdiag;
end


%%%%% recover each fibers using PCA
% coefficients
dt = 1/(M-1);
for ncp = 1:nocomp % # of components used to reconstruct fibers
    
    % for each fiber
    for i=1:N
        cfiber = normalized_f(:,:,i);
        %calculate the pca coefficient
        reconf = zeros(dim,M);
        for pc = 1:dim
            for jj=1:ncp
                tmp_pca_coeff(pc,jj) = dot(squeeze(cfiber(pc,:))',squeeze(pca_basis(pc,:,jj)));
                reconf(pc,:) = reconf(pc,:) + tmp_pca_coeff(pc,jj)*squeeze(pca_basis(pc,:,jj));
            end
        end
        diff(ncp,i) = sqrt(sum(sum((cfiber - reconf).^2,1))*dt);
        
    end
end

%average distance between recovered and original fibers. 
fdiff = mean(diff,2);


