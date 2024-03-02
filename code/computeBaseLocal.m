function [basePixWiseFine, basePixWiseCoarse] = computeBaseLocal(imgCoarse)
% Compute time-series base for each fine pixel
% First on coarse window, then assign to fine pixels
global WINDOW_SIZE_HALF
global N_PIX_COARSE N_TIME
[nlc, nsc, ~] = size(imgCoarse);

basePixWiseCoarse = cell(N_PIX_COARSE, 1);

for ind = 1:N_PIX_COARSE
    [i, j] = ind2sub([nlc,nsc],ind);
    idc = max(1, i - WINDOW_SIZE_HALF) : min(nlc, i + WINDOW_SIZE_HALF);
    jdc = max(1, j - WINDOW_SIZE_HALF) : min(nsc, j + WINDOW_SIZE_HALF);
    samples = reshape(imgCoarse(idc,jdc,:), [], N_TIME)';
    
    basis_i = compute_basis(samples);
    basePixWiseCoarse(ind) = {basis_i};
    
end

% assign to each fine pixel (memory-costly)
basePixWiseFine = assign_basis_c2f(basePixWiseCoarse);

end

function [basis, p_var_cum] = compute_basis_svd(M)
% Input:
%	M: m-by-n, m = length of time series, n = number of samples

[U, var, ~] = svd(M, 'econ');
basis = U * var;
var = diag(var);

sum_var = sum(var);
p_var = var / sum_var;
p_var_cum = cumsum(p_var);

end

function [basis_use] = compute_basis(M)
% Extract basis from samples

global THRES_P_VAR_CUM

[basis, p_var_cum] = compute_basis_svd(M);
n_basis_use = find(p_var_cum >= THRES_P_VAR_CUM, 1);

mag = floor(log10( max(abs(basis(:)))));
basis_use = cat(2, basis(:, 1:n_basis_use), ones(size(basis,1), 1) * 10^mag);

end

function [basis_pix_wise_fine] = assign_basis_c2f(basis_pix_wise_coarse)

global SIZE_COARSE SIZE_FINE N_PIX_FINE SCALE_FACTOR
basis_pix_wise_fine = cell(N_PIX_FINE, 1);

for i = 1:SIZE_COARSE(1)
    rng_if = (i-1)*SCALE_FACTOR +1 : i*SCALE_FACTOR;
    for j = 1:SIZE_COARSE(2)
        rng_jf = (j-1)*SCALE_FACTOR +1 : j*SCALE_FACTOR;
        [X,Y] = meshgrid(rng_if, rng_jf);
        ind_fine_local = reshape(sub2ind(SIZE_FINE,X,Y)', [], 1);    
        ind_coarse = sub2ind([SIZE_COARSE(1),SIZE_COARSE(2)], i, j);
        
        basis_pix_wise_fine(ind_fine_local) = basis_pix_wise_coarse(ind_coarse);

    end
end
    
end
