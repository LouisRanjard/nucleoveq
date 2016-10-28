function [weight_trimmed ] = entrop_trim(weight,threshold)
% trim a weight sequence based on the entropy value and the coverage 
%
% weight: 4 columns sequence of nucleotides
% 

    [ meanEntropy, shannonEnt ] = shannonEntropy(weight) ;
    
    if nargin<2
        threshold=meanEntropy;
    end

    % find the longest contiguous segment
    x = diff([false shannonEnt<threshold false]);
    p = find(x==1);
    q = find(x==-1);
    %[maxlen,ix] = max(q-p);
    [~,ix] = max(q-p);
    first = p(ix);
    last = q(ix)-1;
    %fprintf(1,'\nLargest chunk is %d long ([%d-%d])\n',maxlen, first, last);
    
    weight_trimmed = weight(:,first:last) ;
    
end