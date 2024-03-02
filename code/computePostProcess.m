function [imgFineFusion] = computePostProcess(imgCoarse, imgFine, imgMaskValid, imgFineFusion)
% A simple post-processing. Sophisticated ones like SG filter may be better
global MAX_VALUE
global SCALE_FACTOR
is_anomaly = imgFineFusion < -1*MAX_VALUE | imgFineFusion > MAX_VALUE;
imgCoarseDS = imresize(imgCoarse, SCALE_FACTOR);
imgFineFusion(is_anomaly) = imgCoarseDS(is_anomaly);
imgFineFusion(imgMaskValid) = imgFine(imgMaskValid);
end