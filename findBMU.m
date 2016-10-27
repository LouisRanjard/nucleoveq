function [ treepath, aligned_pos ] = findBMU(seqvect, tree, weight, ce, fe, mode)
% find Best Matching Unit, starting at root of tree and 
% if mode="full" then return the full path with scores
% else only return the BMU tip id and score

    if nargin<6
        mode='';
        if nargin<5
            fe=0;
            if nargin<4
                ce=0;
            end
        end
    end
    
    treepath = zeros(1,2) ;
    treepath(1,:) = [1 nan] ; % starts at the top of the tree, the distance to the root weight is not calculated
    dep = 1 ;

    if (strfind(mode,'full')) % store the full path to BMU
        while (sum(tree==treepath(dep,1))>0) % this node is a father => keep on searching
            leaves = find(tree==treepath(dep,1)) ;
            dep = dep + 1 ;
            scor = zeros(1,max(leaves))*nan ; % initialise to NaN because the distance to the root weight is not calculated and unknown
            for L=leaves
                [ scor(L), ~, aligned_pos ] = DTWaverage( weight{L}, seqvect, 1, 0, ce, fe ) ; % coefficient alignment
            end
            treepath(dep,2) = min(scor(leaves)) ;
            %fprintf(1,'%f\n',BMU(s,dep,2)) ;
            treepath(dep,1) = datasample( find(scor==min(scor(leaves))) ,1 ) ; % find best match (Best Matching Unit), minimum distance
        end
    else % only store last BMU index and score
        while (sum(tree==treepath(1,1))>0) % this node is a father => keep on searching
            leaves = find(tree==treepath(1,1)) ;
            scor = zeros(1,max(leaves))*nan ; % initialise to NaN because the distance to the root weight is not calculated and unknown
            for L=leaves
                [ scor(L), ~, aligned_pos ] = DTWaverage( weight{L}, seqvect, 1, 0, ce, fe ) ; % coefficient alignment
            end
            treepath(1,2) = min(scor(leaves)) ;
            treepath(1,1) = datasample( find(scor==min(scor(leaves))) ,1 ) ; % find best match (Best Matching Unit), minimum distance
        end
    end
    
end
