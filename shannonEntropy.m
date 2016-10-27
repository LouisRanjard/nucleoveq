function [ meanEntropy, shannonEnt ] = shannonEntropy(weightmat)
% compute Shannon entropy for a given weight matrix
% average the entropy value of all positions (columns 1 to 4) in weightmat
%
% For 4 rows: 
% 0 <= shannonEnt <= -sum([.25 .25 .25 .25].*log([.25 .25 .25 .25]))=1.3863
%

    len = size(weightmat,2) ;
    shannonEnt = zeros(1,len) ;
    for n=1:len
        non_zeros_freq = weightmat(weightmat(1:4,n)>0,n) ;
        shannonEnt(n) = -sum(non_zeros_freq.*log(non_zeros_freq)) ;
    end
    
    meanEntropy = sum(shannonEnt)/len ;
    
end
