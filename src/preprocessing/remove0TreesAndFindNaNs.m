function [FVc, nans, species] = remove0TreesAndFindNaNs(FVc, species)
%REMOVE0TREESANDFINDNANS Removes trees that have only 0 values for features
%and finds NaNs from feature data.
%   Removes trees that have a 0 value for all their features (not useful
%   information) and make logical array of places where there are NaNs in
%   the feature data.
nans = zeros(size(FVc,1),size(FVc,2));
stayingTrees = true(1,size(FVc,2));
for i=1:size(FVc,2)
    if FVc(:,i) == 0
        stayingTrees(i) = false;
    else
        nans(:,i) = isnan(FVc(:,i));
    end  
end
FVc = FVc(:,stayingTrees);
nans = nans(:,stayingTrees);
species = species(stayingTrees);
end

