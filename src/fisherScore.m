function [scores] = fisherScore(X,Y,feat,gamma)
%FISHERSCORE Calculates fisher score for each feature of X individually, on
%the condition that features feat have already been included in the reduced
%data space.
%   Rows of X correspond to features, columns to observations. Set feat=[]
%   to calculate fisher scores for single features.
scores = zeros(size(X,1),2);
if isempty(feat)   
    unqSpecies = unique(Y);
    m = arrayfun(@(x) sum(Y == x), unqSpecies);
    for i=1:size(X,1)
        numerator = 0;
        denominator = 0;
        for j=1:size(unqSpecies,2)
            class_j = Y == unqSpecies(j);
            classmean = mean(X(i,class_j),'omitnan');
            numerator = numerator + m(j)*(classmean - mean(X(i,:),'omitnan'))^2;
            denominator = denominator + m(j)*std(X(i,class_j),'omitnan')^2;
        end
        scores(i,:) = [i,numerator/denominator];
    end
else
    X = replaceNaNandInf(X,Y);
    for i=1:size(X,1)
        % If feature i is already included, continue to next iteration
        if ~isempty(find(feat == i,1)) 
            scores(i,:) = [i,0];
            continue;
        end
        tempsubset = [feat, i];
        [B,~,T] = scattermat(X(tempsubset,:),Y);
        scores(i,:) = [i,trace(B/(T + gamma*eye(size(T))))];
    end
end

