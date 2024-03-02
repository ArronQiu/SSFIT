function [imgCoarseRef, imgRMax] = computeCoarseRef(imgCoarse, imgFine, imgMaskValid)
% For each fine pixel, compute the most similar coarse series
% Output
%   imgCoarseRef: the most similar coarse series
%   imgRMax: r of the most similar coarse series
global SIZE_FINE N_PIX_FINE
global THRES_R_CR WINDOW_SIZE_HALF_CR

% Using global variable in PARFOR may cause errors
thres_r_cr = THRES_R_CR;
window_size_half_cr = WINDOW_SIZE_HALF_CR;

imgCoarseRef = zeros(N_PIX_FINE, SIZE_FINE(3), class(imgFine));
imgRMax = zeros(N_PIX_FINE, 1);

for indFine = 1:N_PIX_FINE
    [imgCoarseRef(indFine,:),imgRMax(indFine)] = computeCoarseRef_core(...
        imgCoarse, imgFine, imgMaskValid, indFine, thres_r_cr, window_size_half_cr);
end

imgCoarseRef = reshape(imgCoarseRef, SIZE_FINE);
imgRMax = reshape(imgRMax, SIZE_FINE(1), SIZE_FINE(2));

end

function [coarseRef, rMax] = computeCoarseRef_core(imgCoarse, imgFine, imgMaskValid, indFine, THRES_R_CR, WINDOW_SIZE_HALF_CR)
% For a fine pixel (indFine), compute the most similar coarse series
% Output:
%   coarseRef: the most similar coarse series
%   rMax: r of the most similar coarse series

% Using global variable in PARFOR may cause errors
SIZE_FINE = size(imgFine);
SCALE_FACTOR = size(imgFine,1)/size(imgCoarse,1);

maskValid = reshape(imgMaskValid, [], SIZE_FINE(3));

[iF, jF] = ind2sub([SIZE_FINE(1), SIZE_FINE(2)], indFine);
iC = fix((iF-1)/SCALE_FACTOR) +1;
jC = fix((jF-1)/SCALE_FACTOR) +1;

isValid_i = maskValid(indFine,:);
fine = reshape(imgFine(iF,jF,:), SIZE_FINE(3), 1);
fine(~isValid_i) = nan;

[coarseWin] = get_coarse_win(imgCoarse, iC, jC, WINDOW_SIZE_HALF_CR);

if sum(isValid_i) <= 2 % cannot compute r due to too few fine observation
    % mean of all neighboring coarse pixels
    coarseRef = mean(coarseWin, 2);
    rMax = 1;
    
else
    % compute r for each neighbor coarse series
    coarseRegSample = coarseWin;
    fineRegSample = fine;
    r = corr(fineRegSample, coarseRegSample, 'rows', 'pairwise');
    
    [rMax, indRMax] = max(r);
    if rMax >= THRES_R_CR                       % the most similar
        coarseRef = coarseWin(:, indRMax(1));
    else                                        % mean
        isToAverage = r>0;
        if sum(isToAverage) == 0, isToAverage = true(size(r)); end
        w = r(isToAverage)';
        w = w ./ sum(w);
        coarseRef = coarseWin(:,isToAverage) * w;
        rMax = mean(r(isToAverage));
        
    end
end
end

function [coarseWin] = get_coarse_win(imgCoarse, i, j, w_half)
% coarseWin: (SIZE_FINE(3), n_pix)
SIZE_COARSE = size(imgCoarse);

idc = max(1, i -w_half) : min(SIZE_COARSE(1), i + w_half);
jdc = max(1, j -w_half) : min(SIZE_COARSE(2), j + w_half);

coarseWin = imgCoarse(idc,jdc,:);
coarseWin = reshape(coarseWin, [], SIZE_COARSE(3))';

end
