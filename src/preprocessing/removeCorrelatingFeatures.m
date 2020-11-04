function [FVc,FNc] = removeCorrelatingFeatures(FVc,FNc,corrlimit)
%REMOVECORRELATINGFEATURES Removes strongly correlating features.
%   Removes features that have greater correlation than corrlimit with some
%   other feature.
correlations = corrcoef(FVc');
strongcorrs = abs(correlations) > corrlimit;
strongcorrs = strongcorrs - eye(size(strongcorrs));

% Figure out how many other features does each feature correlate with
% and starting with the one that correlates strongly with most other 
% features (and in descending order after that), remove all features that 
% strongly correlate with it
corrfeatures = sum(strongcorrs,2);
[maxcorrs, i] = max(corrfeatures);
stayingfeatures = true(1,size(FNc,1));
while maxcorrs > 0
    j = 1;
    while j <= size(strongcorrs(i,:),2)
        if strongcorrs(i,j) == 1
            stayingfeatures(j) = false;   
        end
        j = j + 1;
    end
    corrfeatures = sum(strongcorrs.*(stayingfeatures'*stayingfeatures), 2);
    [maxcorrs, i] = max(corrfeatures);
end
FVc = FVc(stayingfeatures,:);
FNc = FNc(stayingfeatures);
end

