function [scores] = constraintScore(X,Y,m,c,lambda)
%CONSTRAINTSCORE Computes the Constraint Score for each feature (row) in X.
%   Y is the responce variable, m is the desired number of must-link pairs
%   for every species, c is the desired number of cannot-link pairs for 
%   every possible pair of different species and lambda is the
%   regularization coefficient used in Constraint Score-2 (Set lambda = -1
%   to compute Constraint Score-1).
unqSpecies = unique(Y);

% Get indices of each species in Y, with same species on the same row
speciesrows = zeros(size(unqSpecies,2),size(Y,2));
for i=1:size(unqSpecies,2)
    speciesrows(i,:) = Y == unqSpecies(i);
end

% Store indices of must-link and cannot-link pairs in Y
mlpairs = zeros(m*size(unqSpecies,2),2);
clpairs = zeros(c*nchoosek(size(unqSpecies,2),2),2);
mlidx = 1;
clidx = 1;
for i=1:size(unqSpecies,2)
    I = find(speciesrows(i,:));
    for j=i:size(unqSpecies,2)
        if i == j
            for k=1:m
                mlpair = randsample(I,2);
                % If the chosen observations have NaN values, choose again
                while nnz(isnan(X(:,mlpair))) ~= 0
                    mlpair = randsample(I,2);
                end
                mlpairs(mlidx,:) = mlpair;
                mlidx = mlidx + 1;
            end
        else
            J = find(speciesrows(j,:));
            for k=1:c
                clpair = [I(randi(numel(I))),J(randi(numel(J)))];
                % If the chosen observations have NaN values, choose again
                while nnz(isnan(X(:,clpair))) ~= 0
                    clpair = [I(randi(numel(I))),J(randi(numel(J)))];
                end
                clpairs(clidx,:) = clpair;
                clidx = clidx + 1;
            end
        end
    end
end

% Compute the sums used to calculate constraint scores
mlsum = zeros(size(X,1),1);
clsum = zeros(size(X,1),1);
for i=1:size(X,1)
    for j=1:size(mlpairs,1)
        mlsum(i) = mlsum(i) + (X(i,mlpairs(j,1)) + X(i,mlpairs(j,2)))^2;
    end
    for j=1:size(clpairs,1)
        clsum(i) = clsum(i) + (X(i,clpairs(j,1)) + X(i,clpairs(j,2)))^2;
    end
end

% Compute constraint scores
if lambda == -1
    scores = [(1:size(X,1))',mlsum./clsum];
else
    scores = [(1:size(X,1))',mlsum - lambda*clsum];
end
end

