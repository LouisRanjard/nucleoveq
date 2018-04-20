function [average_minimum_dist, pct_mini_is_true] = avedist(tree, weight, true_seq)
% Compute the shortest distance similarity between weights and true sequences
% Remove each true sequence when it has been assigned a weight matrix
 
    % if a tree is input then only consider the tips
    if (numel(tree)>1)
        tips = setdiff(1:length(tree),tree);
    else % otherwise consider all weight matrices
        tips = 1:numel(weight) ;
    end
    
    if iscell(weight)
        weight_seq=[];
        s=1;
        for n=tips
            weight_seq(s).Header = ['weight' num2str(n)] ;
            weight_seq(s).Sequence = mat2nucleo(weight{n}) ;
            weight_seq(s).seqvect = weight{n}(1:4,:) ; % get rid of persistence vector
            s=s+1;
        end
    else
        weight_seq = weight;
    end
    
    % make sure the data structure dimension are consistent
    if size(weight_seq,2)>size(weight_seq,1), weight_seq=weight_seq'; end
    if size(true_seq,2)>size(true_seq,1), true_seq=true_seq'; end
    
    % compute percentage similarity between weight and true sequences
    d = zeros(numel(weight_seq),numel(true_seq)) ;
    for n = 1:numel(weight_seq)
        for m = 1:numel(true_seq)
            d(n,m) = seqpdist([weight_seq(n);true_seq(m)],'Method','p-distance') ;
        end
    end
    
    % for each weight find the closest of the true sequences
    minid=ones(1,numel(weight_seq)) ;
    for n = 1:numel(weight_seq)
        [M,I] = min(d(:)) ;
        minid(n) = M ;
        [I_row, I_col] = ind2sub(size(d),I) ;
        d(I_row,:) = 1 ; % delete the weight M distances
        d(:,I_col) = 1 ; % delete the true sequence N as it has already been chosen
    end

    % return the average percentage similarity
    average_minimum_dist = mean(1-minid) ;
    
    % compute percentage similarity between all sequences and check if weight are pulled toward a true_seq (and not another weight)
    % [duplication: could be done at the same time as above rather than recomputing the distances]
    da = seqpdist([weight_seq;true_seq],'Method','p-distance') ;
    das = squareform(da) ;
    das(logical(eye(size(das)))) = +Inf ; % set the diagonal (identity) to +Infinity instead of 0
    mini_is_true = zeros(1,numel(weight_seq)) ;
    for n = 1:numel(weight_seq)
        [~,I] = min(das(n,:)) ; % find the closest sequence to weight n
        if I>n % is it another weight or a true_seq
            mini_is_true(n)=1 ;
        end
    end
    
    % return the percentage of weight sequence which have a true sequence as their closest 
    pct_mini_is_true = sum(mini_is_true)/numel(weight_seq) ;
    
end
