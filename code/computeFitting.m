function [imgFineFusion] = computeFitting(fineFill, basePixWiseFine, weight)
% Fit the filled series using extracted bases and the weight
global SIZE_FINE
global CLASS_FINE
fineFit = zeros(size(fineFill),CLASS_FINE);
n_to_fit = size(fineFill, 1);
for i = 1:n_to_fit
    y = single(fineFill(i, :)');
    A = single(basePixWiseFine{i});
    W = single(diag(weight(i, :)));
    AW = W*A;
    yW = W*y;
    yfit = A * (pinv(AW)*yW);
    fineFit(i, :) = cast(yfit, CLASS_FINE);
    
end
imgFineFusion = reshape(fineFit, SIZE_FINE);
end

