function [subset,bestloss] = addToImproveAccuracy(FV,species,subset,bestloss,featureOrder,ncand,minimpr,maxtime,NameValueArgs)
%ADDTOIMPROVEACCURACY Add a feature to a subset if it improves
%classification accuracy
%   Add features to subset one by one according to featureOrder to see if 
%   classification accuracy improves. Finds ncand number of candidate
%   subsets that have at least minimpr lower loss than bestavgloss and 
%   returns the one with the lowest loss. Specify maxtime (in seconds) to 
%   terminate execution after a certain amount of time has passed after the
%   last improvement.
arguments
    % Feature values
    FV
    species
    subset
    bestloss
    % Ranking of features acquired from a filter
    featureOrder
    % Parameters of DIE1VA
    ncand
    minimpr
    maxtime
    NameValueArgs.CMethod string = 'knn' 
    NameValueArgs.NumNeighbors uint32 = 1
    NameValueArgs.Distance string = 'euclidean'
    NameValueArgs.Weighed string = 'off'
end

% Temporary subset where we add new features to the already selected feature set
tempsubset = [subset,0];

% Matrix and vector for saving candidate subsets and their losses
candidatesubsets = zeros(ncand,size(tempsubset,2));
betterlosses = ones(ncand,1);

% The index where we save new candidate subsets and their losses
k = 1;

% Apply balanced or default class weights
if strcmp(NameValueArgs.Weighed,'on')
    W = applyWeights(species);
else
    W = ones(length(species),1);
end

% Start the clock and go through features in the order of their ranking
elapsedTime = 0;
for i=1:size(featureOrder,2)
    tic
    
    % If feature has not been selected yet, add it to tempsubset
    if ~ismember(featureOrder(i),tempsubset)
        tempsubset(end) = featureOrder(i);
    
    % If feature has already been selected, move on to the next feature
    else
        continue;
    end

    % Build classifier
    if strcmp(NameValueArgs.CMethod,'knn')
        mdl = fitcknn(FV(tempsubset,:)',species,'KFold',10, ...
            'NumNeighbors',NameValueArgs.NumNeighbors,'Weights',W, ...
            'Distance',NameValueArgs.Distance);
    elseif strcmp(NameValueArgs.CMethod,'svm')
        mdl = fitcsvm(FV(tempsubset,:)',species,'KFold',10,'Weights',W);
    end

    % Evaluate loss
    loss = kfoldLoss(mdl);
    
    % If improvement in accuracy is more than minimpr, save new candidate subset
    if loss < bestloss - minimpr

        % Reset the clock, since an improvement was found
        elapsedTime = 0;

        % Save new candidate subset and corresponding loss
        betterlosses(k) = loss;
        candidatesubsets(k,:) = tempsubset;

        % If there are ncand candidate subsets, select the one with the lowest
        % loss and stop searching
        if k == ncand
            [bestloss,besti] = min(betterlosses);
            subset = candidatesubsets(besti,:);
            break;
        end

        % Move on to next index
        k = k + 1;
    end

    % Keep track of elapsed time and stop searching if it exceeds maxtime
    elapsedTime = elapsedTime + toc;
    if maxtime < elapsedTime
        break;
    end
end

% If any candidate subsets were found, choose the one with the lowest loss
% (betterloss == 1 if no improvements were found)
[betterloss,besti] = min(betterlosses);
if betterloss ~= 1
    bestloss = betterloss;
    subset = candidatesubsets(besti,:);
end
end