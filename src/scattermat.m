function [B, W, T] = scattermat(X,Y)
%SCATTERMAT Calculate scatter matrices (between-class, within-class and
%total)
%   B:Between-class scatter matrix
%   W:Within-class scatter matrix
%   T:Total scatter matrix
%   data: data space spanned by selected features (features as rows,
%   columns as observations)
%   Y: class information of observations
%
    [l, ~] = size(X);
    classes = unique(Y);
    tot_classes = length(classes);
    B = zeros(l,l);
    W = zeros(l,l);           
    overallmean = mean(X,2);
    for i=1:tot_classes
        classi = find(Y == classes(i));
        xi = X(:,classi)';
        
        mci = mean(xi);
        xi = xi-repmat(mci,length(classi),1);
        W = W + xi'*xi;
        B = B + length(classi)*(mci-overallmean)*(mci-overallmean)';
    end
    T = zeros(l,l);
    for i=1:length(Y)
        zi = X(:,i);
        T = T + (zi - overallmean)*(zi - overallmean)';
    end
end

