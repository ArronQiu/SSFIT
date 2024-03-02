function [fineFill, weight] = computeFillSimple(imgCoarse, imgFine, imgMaskValid)
% To ensure enought fitting points and robustness, particularly when there
% are insufficient fine-resolution observations, gap-fill before fitting.
% Note that, the filling here is quite simple, sophisticated ones like NSPI
% may be worthy of exploration.

% step 1: find coarse reference series
[imgCoarseRef, imgRMax] = computeCoarseRef(imgCoarse, imgFine, imgMaskValid);
global N_TIME
% step 2: fill by correcting systematic bias
fine = reshape(imgFine, [], N_TIME);
coarseRef = reshape(imgCoarseRef, [], N_TIME);
maskValid = reshape(imgMaskValid, [], N_TIME);
fineFill = computeFill(fine, coarseRef, maskValid);

% also, compute the weight of the filling points in fitting
global N_PIX_FINE
rMax = reshape(imgRMax, N_PIX_FINE, 1);
weight = computeWeight(rMax, maskValid);

end

function [fineFill] = computeFill(fine, coarseRef, maskValid)
bias = computeBiasSystematic(fine, coarseRef, maskValid);
fineFill = coarseRef + bias;

% avoid anomaly
isOutBound = fineFill > 1e4 | fineFill< -1e4;
fineFill(isOutBound) = coarseRef(isOutBound);

% keep original observations
fineFill(maskValid) = fine(maskValid);

end

function [bias] = computeBiasSystematic(fine, coarseDS, maskValid)
% compute systematic bias based on valid observations at fine resolution
fine(~maskValid) = nan;
coarseDS(~maskValid) = nan;

fine_mean = mean(fine, 2, 'omitnan');
coarseDS_mean = mean(coarseDS, 2, 'omitnan');

bias = fine_mean - coarseDS_mean;

end

function  [weight] = computeWeight(rMax, maskValid)
% compute the weight in fitting
global N_PIX_FINE SIZE_FINE N_TIME
weight = ones(N_PIX_FINE, SIZE_FINE(3));
rMax = repmat(rMax, 1, N_TIME);
weight(~maskValid) = rMax(~maskValid);
end
