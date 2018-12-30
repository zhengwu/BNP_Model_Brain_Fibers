function final_reconstl = roifibers_reconstruction_adaptive(labels,Trans,Rots,PCA_coeffs,means_sl,pca_basis)
% function for reconstruction of fibers connecting two rois.
Nfibs = length(labels);
dim = size(means_sl,1);
%[dim,nocomp,Nfibs] = size(PCA_coeffs);
[~,~,M,~]=size(pca_basis);

%reconstruct the signal;
recon_f = zeros(dim,M,Nfibs);
for i=1:Nfibs
    label = labels(i);
    currt_mn = means_sl(:,:,label);
    pca_coeff = PCA_coeffs{i};
    for pc=1:dim
        tmp_pca_coeffs = pca_coeff{pc};
        nocomp = length(tmp_pca_coeffs);
        for jj=1:nocomp
            recon_f(pc,:,i) = recon_f(pc,:,i) + tmp_pca_coeffs(jj)*squeeze(pca_basis(label,pc,:,jj))';
        end;
        recon_f(pc,:,i) = recon_f(pc,:,i) + currt_mn(pc,:);
    end;
end;

%recover each with translation and rotation
final_reconstl = zeros(dim,M,Nfibs);
for i=1:Nfibs
    final_reconstl(:,:,i) = squeeze(Rots(:,:,i))'*squeeze(recon_f(:,:,i)) + squeeze(Trans(:,i))*ones(1,M);
end;
