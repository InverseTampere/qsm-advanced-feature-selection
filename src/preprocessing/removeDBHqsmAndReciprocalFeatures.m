function [FV,FN] = removeDBHqsmAndReciprocalFeatures(FV,FN)
%REMOVEDBHQSMANDRECIPROCALFEATURES Removes features that have 'DBHqsm' in
%them and reciprocals of other features.
%   Removes features with DBHqsm in them (as it is almost always exactly 
%   the same as DBHcyl), as well as reciprocals of other features which
%   don't correlate linearly with other features but are dependent on them
%   and as such, useless for the classifier.
k = strfind(FN,'DBHqsm');
keepfeatures = true(size(FN,1),1);
for i=1:size(FN,1)
    if ~isempty(k{i})
        keepfeatures(i) = false; 
    end
end
FN = FN(keepfeatures);
FV = FV(keepfeatures,:);
keepfeatures = true(size(FN,1),1);
for i=16:size(FN,1)
    [numerator, denominator] = strtok(FN{i},'/');
    if ~isempty(denominator)
        denominator = denominator(2:end);
        reciprocal = denominator + "/" + numerator;
        for j=i:size(FN,1)
            if string(FN{j}) == reciprocal
                keepfeatures(j) = false;
                break;
            end
        end
    end    
end
FN = FN(keepfeatures);
FV = FV(keepfeatures,:);
end

