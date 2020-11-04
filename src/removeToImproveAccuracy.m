function [subset,loss] = removeToImproveAccuracy(FV,species,subset,loss,NameValueArgs)
%REMOVETOIMPROVEACCURACY Remove a feature if it improves accuracy.
%   Removes the one feature that improves classification accuracy the most
%   (if any). Rows of FV correspond to features, columns
%   to observations.
arguments
    % Feature values
    FV
    species
    subset
    loss
    NameValueArgs.CMethod string = 'knn'
    NameValueArgs.NumNeighbors uint32 = 1
    NameValueArgs.Distance string = 'euclidean'
    NameValueArgs.Weighed string = 'off'
end

len = length(subset);

% Create all possible feature combinations that can be made by removing one
% feature from the subset, initialize vector to save the corresponding losses
featuresubsets = nchoosek(subset,len-1);
iter = size(featuresubsets,1);
losses = zeros(iter,1);

% Loop through every possible combination
for k=1:iter
    
    inputspace = FV(featuresubsets(k,:),:);

    % Create partition and initialize the number of incorrect predictions
    cv = cvpartition(length(species),'KFold',10); 
    wrong = 0;

    % Loop through test sets in evaluating the classifier
    for i=1:cv.NumTestSets

        % Get training and test sets
        trainSet = cv.training(i);
        testSetIdx = find(cv.test(i));

        % Apply balanced or default class weights on the training set
        if strcmp(NameValueArgs.Weighed,'on')
            W = applyWeights(species(trainSet));
        else
            W = ones(nnz(trainSet),1);
        end

        % Build classifier based on training set
        if strcmp(NameValueArgs.CMethod,'knn')
            mdl = fitcknn(inputspace(:,trainSet)',species(trainSet), ...
                'NumNeighbors',NameValueArgs.NumNeighbors,'Weights',W, ...
                'Distance',NameValueArgs.Distance);
        elseif strcmp(NameValueArgs.CMethod,'svm')
            mdl = fitcsvm(inputspace(:,trainSet)',species(trainSet),...
                'Weights',W);
        end

        % Predict the class of each observation in test set
        for j=1:cv.TestSize(i)
            label = predict(mdl,inputspace(:,testSetIdx(j))');

            % If prediction was incorrect, update the number of incorrect
            % predictions
            if label ~= species(testSetIdx(j))
                wrong = wrong + 1;
            end
        end
    end

    % Save the loss of this feature set
    losses(k) = wrong/length(species);
end

% Get the subset with the lowest loss
[newloss,besti] = min(losses);
newsubset = featuresubsets(besti,:);

% If new loss is less than the old one, update subset and loss
if newloss < loss
    subset = newsubset;
    loss = newloss;
end
end

