function [ meanEntropy, shannonEnt, averageent ] = shannonEntropy_s(weightmats,normalise,averageent)
% compute Shannon entropy for the consensus sequence of a set of weight matrices
%

    if nargin<3
        averageent=0;
        if nargin<2
            normalise=1;
        end
    end
    
    % compute the average entropy of each weight matrix independently
    if averageent>0
        averageent=mean(cellfun(@shannonEntropy,weightmats,'UniformOutput',0));
    end
    
    % normalise: convert each weightmat to a minimum entropy sequence
    % each position takes only 1 and 0 according to the frequency of each base (same way weight sequences are defined)
    if normalise
        weightmatsn = cellfun(@round,weightmats,'UniformOutput',0) ;
    else
        weightmatsn = weightmats ;
    end
    
    summat = {zeros(size(weightmatsn{1},1),size(weightmatsn{1},2))} ;
    
    % element wise sum
    for n=1:numel(weightmatsn)
        summat{1}(1:4,:) = summat{1}(1:4,:) + weightmatsn{n}(1:4,:) ;
    end
    
    % normalise by number of weight to get frequencies
    summat{1}(1:4,:) = summat{1}(1:4,:)./numel(weightmatsn) ;
    
    [ meanEntropy, shannonEnt ] = shannonEntropy(summat{1}) ;
    
end