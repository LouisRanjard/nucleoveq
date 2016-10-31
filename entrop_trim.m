function [weight_trimmed ] = entrop_trim(weight,threshold,momean)
% trim a weight sequence based on the entropy value and the coverage 
%
% weight: 4 columns sequence of nucleotides
% 

    if nargin<3
        momean=0 ;
    end

    if iscell(weight)
        %meanEntropy = cellfun(@(x) shannonEntropy(x), weight) ;
        meanEntropy = zeros(numel(weight),1) ;
        shannonEnt = cell(numel(weight),1) ;
        for w = 1:numel(weight)
            [ meanEntropy(w), shannonEnt{w} ] = shannonEntropy(weight{w}) ;
        end
        meanEntropy = mean(meanEntropy) ;
    else
        [ meanEntropy, shannonEnt ] = shannonEntropy(weight) ;
    end
    
    if nargin<2
        threshold = meanEntropy ;
    end
    
    if iscell(weight)
        weight_trimmed = cell(numel(weight),1) ;
        for w = 1:numel(weight)
            if momean
                [first,last] = ent_trim(movimean(shannonEnt{w}),threshold) ;
            else
                [first,last] = ent_trim(shannonEnt{w},threshold) ;
            end
            weight_trimmed{w} = weight{w}(:,first:last) ;
        end
    else
        if momean
            [first,last] = ent_trim(movimean(shannonEnt),threshold) ;
        else
            [first,last] = ent_trim(shannonEnt,threshold) ;
        end
        weight_trimmed = weight(:,first:last) ;
    end
    
end

%% local function
function [ first, last ] = ent_trim(entrop,thresh)
    % find the longest contiguous segment of entropy that is greater than threshold
    x = diff([false entrop<=thresh false]);
    p = find(x==1);
    q = find(x==-1);
    %[maxlen,ix] = max(q-p);
    [~,ix] = max(q-p);
    first = p(ix);
    last = q(ix)-1;
    %fprintf(1,'\nLargest chunk is %d long ([%d-%d])\n',maxlen, first, last);
end
