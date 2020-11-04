function [FVc] = replaceNaNandInf(FVc,species)
%REPLACENANANDINF Replaces NaN and Inf values in matrix FVc.
%   Replaces NaN and Inf values with the mean and finite maximum values
%   of that particular feature for that particular species, respectively.
for i=1:size(FVc,1)
    for j=1:size(FVc,2)
        if isnan(FVc(i,j))
            speciestrees = species == species(j);
            FVc(i,j) = mean(FVc(i,speciestrees),'omitnan');
        end
        if isinf(FVc(i,j))
            speciestrees = species == species(j);
            FVc(i,j) = max(isfinite(FVc(i,speciestrees)),[],'omitnan');
        end
    end
end
end

