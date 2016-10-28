function [weight_trimmed ] = entrop_trim(weight)
% trim a weight sequence based on the entropy value and the coverage 
%
% weight: 4 columns sequence of nucleotides
% 

    weight_trimmed = weight ;

    t = [0 2 2 3 0 2 3 0 2 2 2 0 2 5 3] ;

    % find the longest contiguous segment
    x = diff([false t>1 false]);
    p = find(x==1);
    q = find(x==-1);
    [maxlen,ix] = max(q-p);
    first = p(ix);
    last = q(ix)-1;
    fprintf(1,'\nLargest chunk is %d long ([%d-%d])\n',maxlen, first, last);
    
end