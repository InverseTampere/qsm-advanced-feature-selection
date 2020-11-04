function [W] = applyWeights(species)
%APPLYWEIGHTS Apply weights to the species according to their frequencies
%in the sample.
%   Applies weights to each unique species in a sample so that
%   misclassification of each species is inversely proportional to its
%   size in the sample to prevent overemphasizing classification of most
%   numerous species.
unqSpecies = unique(species);
m = arrayfun(@(x) sum(species == x), unqSpecies);
W = zeros(1,size(species,2));
for i=1:size(species,2)
    for j=1:size(unqSpecies,2)
        if species(i) == unqSpecies(j)
            W(i) = 1/m(j);
        end
    end
end
end

