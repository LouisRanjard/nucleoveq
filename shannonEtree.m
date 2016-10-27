function [ wshannon ] = shannonEtree(BMU,tree,weight)
% for each tip in the tree compute the shannon entropy weighted by the
% number of vector matching the node


    tips = setdiff(1:numel(tree),unique(tree)) ;
    wshannon = zeros(1,length(tips)) ;
    
    for t=1:length(tips)
        matches = sum(BMU(:,1)==tips(t)) ;
        match_siblings = sum(BMU(:,1)==find(tree==tree(tips(t)))) ;
        wshannon(t) = shannonEntropy() ;
    end

end
