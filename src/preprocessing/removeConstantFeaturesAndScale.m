function [FVc, FNc] = removeConstantFeaturesAndScale(FVc, FNc)
%REMOVECONSTANTFEATURESANDSCALE Scales features to [0 1] and removes
%features that have constant values.
%   Scales features to [0 1] and removes features that have min == max (for
%   example TrunkVolume/TrunkVolume).
i = 1;
stayingFeatures = true(size(FVc, 1), 1);
while i <= size(FVc,1)
    feature = FVc(i,:);
    fmin = min(feature(isfinite(feature)));
    fmax = max(feature(isfinite(feature)));
    if fmin < fmax
        FVc(i,:) = (FVc(i,:) - fmin)/(fmax-fmin);
    else
        stayingFeatures(i) = false;
    end
    i = i + 1;
end
FVc = FVc(stayingFeatures,:);
FNc = FNc(stayingFeatures,:);
end

