function [ shannonEnt ] = shannonEntropy_s(weightmats)
% compute Shannon entropy for the consensus sequence of a set of weight matrices
%

    % convert each weightmat to a minimum entropy sequence
    weightmatsn = cellfun(@round,weightmats,'UniformOutput',0) ;
    
    summat = {zeros(size(weightmatsn{1},1),size(weightmatsn{1},2))} ;
    
    % element wise sum
    for n=1:numel(weightmatsn)
        summat{1}(1:4,:) = summat{1}(1:4,:) + weightmatsn{n}(1:4,:) ;
    end
    
    % normalise by number of weight to get frequencies
    summat{1}(1:4,:) = summat{1}(1:4,:)./numel(weightmatsn) ;
    
    shannonEnt = shannonEntropy(summat{1}) ;
    
end