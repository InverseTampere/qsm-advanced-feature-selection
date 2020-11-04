function [pred,preds,acc,elapsedTime,finalsubsets,finallosses] = ...
    DIE1VA(X,Y,n,n_HFS,m,n_b,minimpr,n_c,maxtime,NameValueArgs)

arguments
    % Feature values
    X
    % Class labels
    Y
    % Parameters of DIE1VA
    n
    n_HFS
    m
    n_b
    minimpr
    n_c
    maxtime
    NameValueArgs.FSMethod string = 'fisher'
    NameValueArgs.CMethod string = 'knn'
    NameValueArgs.NumNeighbors uint32 = 1
    NameValueArgs.Distance string = 'euclidean'
    NameValueArgs.Weighed string = 'off'
end

tStart = tic;
unqspecies = unique(Y);
finalsubsets = zeros(length(unqspecies)*n,20);
finallosses = zeros(length(unqspecies)*n,1);
for i=1:length(unqspecies)
    
    % Combine other species into one class, 'other'
    otheri = 1:length(unqspecies);
    otheri(i) = [];
    other = categories(removecats(unqspecies(otheri)));
    onevsall = mergecats(Y,other,'other');
    
    % Apply balanced or default class weights
    if strcmp(NameValueArgs.Weighed,'on')
        W = applyWeights(onevsall);
    else
        W = ones(length(Y),1);
    end
    
    % Rank features using the chosen filter
    if strcmp(NameValueArgs.FSMethod,'fisher')
        scores = fisherScore(X,onevsall,[],0);
        scores = sortrows(scores,2,'descend');
        idx = scores(:,1)';
    elseif strcmp(NameValueArgs.FSMethod,'chi2')
        [idx,~] = fscchi2(X',onevsall);
    elseif strcmp(NameValueArgs.FSMethod,'mrmr')
        [idx,~] = fscmrmr(X',onevsall);
    elseif strcmp(NameValueArgs.FSMethod,'relief')
        [idx,~] = relieff(X',onevsall,1);
    elseif strcmp(NameValueArgs.FSMethod,'constraint')
        scores = constraintScore(X,onevsall,50,10,0.1);
        scores = sortrows(scores,2,'descend');
        idx = scores(:,1)';
    end
    
    % Generate feature subsets using HFS and save their losses
    subsets = zeros(n_b,20);
    losses = zeros(n_b,1);
    for j=1:n_b
        [subset,loss] = HybridFS(X,onevsall,idx',n_HFS,m,W, ...
            'CMethod',NameValueArgs.CMethod, ...
            'NumNeighbors',NameValueArgs.NumNeighbors, ...
            'Distance',NameValueArgs.Distance);

        % Check if accuracy can be improved by removing features
        while length(subset) > 1
            [subset,betterloss] = removeToImproveAccuracy(X,onevsall, ...
                subset,loss,'CMethod',NameValueArgs.CMethod, ...
                'Weighed',NameValueArgs.Weighed, ...
                'Distance',NameValueArgs.Distance, ...
                'NumNeighbors',NameValueArgs.NumNeighbors);

            % Move on if accuracy did not improve
            if loss == betterloss
                break;
            end

            % Save new loss if accuracy improved
            loss = betterloss;
        end
        subsets(j,1:length(subset)) = subset;
        losses(j) = loss;
        idx(1:m) = [];
    end
    
    % Sort HFS-generated subsets according to their losses
    subsets = [subsets, losses];
    subsets = sortrows(subsets,size(subsets,2));
    losses = subsets(:,end);
    subsets = subsets(:,1:end-1);
    
    % Refine HFS-generated subsets starting with the one with the lowest
    % loss
    for x=1:n
        subset = subsets(1,subsets(1,:) > 0);
        prevloss = losses(1);
        subsets(1,:) = [];
        
        % Combine HFS-generated subsets to see if accuracy improves enough
        for j=1:size(subsets,1)
            newsubset = union(subset,subsets(j,subsets(j,:) > 0));
            if newsubset(1) == 0
                newsubset(1) = [];
            end
            
            % Build classifier
            if strcmp(NameValueArgs.CMethod,'knn')
                mdl = fitcknn(X(newsubset,:)',onevsall,'KFold',10, ...
                    'Weights',W,'NumNeighbors',NameValueArgs.NumNeighbors, ...
                    'Distance',NameValueArgs.Distance);
            elseif strcmp(NameValueArgs.CMethod,'svm')
                mdl = fitcsvm(X(newsubset,:)',onevsall,'Weights',W, ...
                    'KFold',10);
            end
            
            % Evaluate loss
            loss = kfoldLoss(mdl);

            % Check if accuracy can be improved by removing features
            while length(newsubset) > 1
                [newsubset,betterloss] = removeToImproveAccuracy(X, ...
                    onevsall,newsubset,loss, ...
                    'CMethod',NameValueArgs.CMethod, ...
                    'Weighed',NameValueArgs.Weighed, ...
                    'Distance',NameValueArgs.Distance, ...
                    'NumNeighbors',NameValueArgs.NumNeighbors);

                % Move on if accuracy did not improve
                if loss == betterloss
                    break;
                end

                % Save new loss if accuracy improved
                loss = betterloss;
            end
            
            % If improvement is not at least minimpr for each new feature,
            % it is not considered enough
            if loss < prevloss - minimpr*(length(newsubset) - ...
                    length(subset))
                subset = newsubset;
                prevloss = loss;
            end
        end
        
        % Add or remove features one by one until there's little to no improvement 
        while length(subset) > 1
            
            % remove features if it improves performance
            [subset,loss] = removeToImproveAccuracy(X,onevsall,subset, ...
                prevloss,'Weighed',NameValueArgs.Weighed, ...
                'CMethod',NameValueArgs.CMethod, ...
                'NumNeighbors',NameValueArgs.NumNeighbors, ...
                'Distance',NameValueArgs.Distance);

            % If accuracy did not improve by removing features, try adding them
            if loss == prevloss
                
                % Re-rank features based on their fisher scores with already 
                % selected features
                fisherscores = fisherScore(X,onevsall,subset,0);
                fisherscores = sortrows(fisherscores,2,'descend');
                idx = fisherscores(:,1);
                
                % add a feature if it improves performance
                [subset,loss] = addToImproveAccuracy(X,onevsall, ...
                    subset,prevloss,idx',n_c,minimpr,maxtime, ...
                    'CMethod',NameValueArgs.CMethod, ...
                    'Weighed',NameValueArgs.Weighed, ...
                    'Distance',NameValueArgs.Distance, ...
                    'NumNeighbors',NameValueArgs.NumNeighbors);

                % Move on if accuracy did not improve
                if loss == prevloss
                    break;
                end     
            end

            % Save new loss if accuracy improved
            prevloss = loss;
        end

        % Save the final subset and corresponding loss
        finalsubsets((i-1)*n+x,1:length(subset)) = subset;
        finallosses((i-1)*n+x) = loss;
    end
end

% Make One-vs-All predictions for every species using features obtained for
% each species
preds = categorical(zeros(length(unqspecies)*n,length(Y)));
for i=1:length(unqspecies)
    for x=1:n

        % Combine other species into one class, 'other'
        otheri = 1:length(unqspecies);
        otheri(i) = [];
        other = categories(removecats(unqspecies(otheri)));
        onevsall = mergecats(Y,other,'other');

        % Apply balanced or default class weights
        if strcmp(NameValueArgs.Weighed,'on')
            W = applyWeights(onevsall);
        else
            W = ones(length(Y),1);
        end

        % Build classifier
        if strcmp(NameValueArgs.CMethod,'knn')
            mdl = fitcknn(X(finalsubsets((i-1)*n+x, ...
                finalsubsets((i-1)*n+x,:) > 0),:)', ...
                onevsall,'Weights',W,'KFold',10, ...
                'Distance',NameValueArgs.Distance, ...
                'NumNeighbors',NameValueArgs.NumNeighbors);
        elseif strcmp(NameValueArgs.CMethod,'svm')
            mdl = fitcsvm(X(finalsubsets((i-1)*n+x, ...
                finalsubsets((i-1)*n+x,:) > 0),:)', ...
                onevsall,'Weights',W,'KFold',10);
        end

        % Make One-vs-All predictions
        preds((i-1)*n+x,:) = kfoldPredict(mdl);
    end
end

% Combine One-vs-All predictions to make final prediction based on their
% votes
pred = categorical(zeros(1,length(Y)));
for i=1:size(preds,2)

    % Get all votes that are not 'other' from one observation
    notother = preds(preds(:,i) ~= 'other',i);

    % If only one vote was not 'other', assing final prediction as that class
    if length(notother) == 1
        pred(i) = notother;

    % If all votes were 'other' assign final prediction randomly from the
    % included species
    elseif isempty(notother)
        pred(i) = unqspecies(randi(length(unqspecies)));

    % If there were more than one vote for specific species, assign final
    % prediction as the species with most votes
    else
        [~,~,C] = mode(notother);

        % If there is a tie, assign final prediction randomly from the
        % species with most votes
        pred(i) = C{1}(randi(length(C)));
    end
end

% If there are extra categories that were not in the predictions, remove them
pred = removecats(pred);

% Calculate accuracy of predictions
acc = nnz(pred == Y)/length(Y);

% Record time consumption
elapsedTime = toc(tStart);
end


