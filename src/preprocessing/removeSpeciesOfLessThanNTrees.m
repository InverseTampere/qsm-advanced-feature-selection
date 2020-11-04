function [FVc, species] = removeSpeciesOfLessThanNTrees(FV, QSMs, n)
%REMOVESPECIESOFLESSTHANNTREES Removes species of less than n trees.
%   Removes species that have less than a required amount of trees for the
%   classifier.
FVc = FV;
% Get each tree's species into categorical array
species = categorical(1,size(QSMs,2));
for i=1:size(QSMs,2)
    species(i) = QSMs(i).species;
end

% Get unique species and the amount of trees belonging to each one
unqSpecies = unique(species);
m = arrayfun(@(x) sum(species == x), unqSpecies);

% Keep trees belonging to species that have a reasonable number (10-15?) of 
% trees for classification
singleSpeciesTreesNeeded = n;
treesKept = zeros(1,size(QSMs,2));
for i=1:size(m,2)
    if m(i) >= singleSpeciesTreesNeeded
        treesKept = treesKept + (species == unqSpecies(i));
    end
end
FVc = FVc(:,logical(treesKept));
species = species(logical(treesKept));
end

