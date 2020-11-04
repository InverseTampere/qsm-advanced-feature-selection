function [FVc, nans, species] = removeTreesWithTooManyNaNs(FVc, nans, NaNpercentagelimit, species)
%REMOVETREESWITHTOOMANYNANS Removes trees that have more than
%NaNpercetagelimit NaNs.
%   Removes trees that have NaN-values for too many features to be useful
%   for the classifier.
treeNaNlimit = size(FVc,1)*NaNpercentagelimit;
treeNaNs = sum(nans,1);
i = 1;
stayingTrees = true(1, size(treeNaNs, 2));
while i <= size(treeNaNs, 2)
    if treeNaNs(i) > treeNaNlimit
        stayingTrees(i) = false;
    end
    i = i + 1;
end
FVc = FVc(:,stayingTrees);
nans = nans(:,stayingTrees);
species = species(stayingTrees);
end

