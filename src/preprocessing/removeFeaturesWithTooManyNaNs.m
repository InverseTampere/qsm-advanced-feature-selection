function [FVc, FNc, nans] = removeFeaturesWithTooManyNaNs(FVc, FNc, nans, NaNpercentagelimit)
%REMOVEFEATURESWITHTOOMANYNANS Removes features that have more than
%NaNpercentagelimit NaNs.
%   Removes features that have NaN-values for too many trees to be useful
%   for the classifier.
featureNaNlimit = size(FVc,2)*NaNpercentagelimit;
featureNaNs = sum(nans,2);
i = 1;
stayingFeatures = true(size(featureNaNs, 1), 1);
while i <= size(featureNaNs, 1)
    if featureNaNs(i) > featureNaNlimit
        stayingFeatures(i) = false;
    end
    i = i + 1;
end
FVc = FVc(stayingFeatures,:);
nans = nans(stayingFeatures,:);
FNc = FNc(stayingFeatures);
end

