function [bestsubset,limit] = HybridFS(X,Y,ranking,n_HFS,m,W,NameValueArgs)
%HYBRID2 Filter-wrapper hybrid method of generating feature subsets for 
%classification.
%   Combines feature subsets according to a ranking to improve the accuracy
%   of the resulting feature set.
arguments
    % Feature values
    X
    % class labels
    Y
    ranking
    n_HFS
    m
    % class weights
    W
    NameValueArgs.random = false
    NameValueArgs.CMethod string = 'knn'
    NameValueArgs.NumNeighbors uint32 = 1
    NameValueArgs.Distance string = 'euclidean'
end

% Randomize ranking
if NameValueArgs.random == true
    ranking(1:m,:) = ranking(randperm(m),:);
end

% Build classifier using first feature in ranking
if strcmp(NameValueArgs.CMethod,'knn')
    mdl = fitcknn(X(ranking(1,:),:)',Y,'KFold',10,'NumNeighbors', ...
        NameValueArgs.NumNeighbors,'Weights',W,'Distance',NameValueArgs.Distance);
elseif strcmp(NameValueArgs.CMethod,'svm')
    mdl = fitcsvm(X(ranking(1,:),:)',Y,'KFold',10,'Weights',W);
end

% Set limit for classification accuracy
limit = kfoldLoss(mdl);
bestsubset = ranking(1,:);
evaluatedSubsets = zeros(n_HFS*(2*m-n_HFS),20);
evaluatedSubsets(1,1:length(bestsubset)) = bestsubset;
evali = 2;
d = 2;
while true
    idx = 1;
    subsets = zeros(n_HFS*(m-n_HFS+m-1)/2,d+1);
    for i=1:n_HFS
        for j=i+1:m
            subset = union(ranking(i,:), ranking(j,:));
            subset = subset(subset > 0);
            % Evaluate each combination of features only once
            if length(subset) > size(evaluatedSubsets,2) || ...
                ismember(subset,evaluatedSubsets(:,1:length(subset)),'rows') ~= 0
                continue;
            end
            
            % Build classifier using new combination of features
            if strcmp(NameValueArgs.CMethod,'knn')
                mdl = fitcknn(X(subset(subset > 0),:)',Y,'KFold',10, ...
                    'NumNeighbors',NameValueArgs.NumNeighbors,'Weights',W, ...
                    'Distance',NameValueArgs.Distance);
            elseif strcmp(NameValueArgs.CMethod,'svm')
                mdl = fitcsvm(X(subset(subset > 0),:)',Y,'KFold',10, ...
                    'Weights',W);
            end
            
            % Evaluate accuracy and save the subset if its accuracy is 
            % higher than the limit
            loss = kfoldLoss(mdl);
            evaluatedSubsets(evali,1:length(subset)) = subset;
            evali = evali + 1;
            if loss < limit
                subsets(idx,1:length(subset)) = subset;
                subsets(idx,d+1) = loss;
                idx = idx + 1;
            end       
        end
    end
    % If accuracy did not improve, end algorithm
    if subsets == 0
        break;
    end
    ranking = sortrows(subsets,d+1);
    % Delete possible rows of zeros
    ranking = ranking(any(ranking,2),:);
    if isempty(ranking)
        break;
    end
    if m > size(ranking,1)
        m = size(ranking,1);
    end
    if n_HFS > size(ranking,1)
        n_HFS = size(ranking,1);
    end
    limit = ranking(1,end);
    bestsubset = ranking(1,ranking(1,:) >= 1);
    ranking = ranking(:,1:end-1);
    d = d*2;
end
end

