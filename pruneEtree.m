function [ pruned_tree, pruned_weight ] = pruneEtree(reads, tree, weight, ce, fe)
% find the tip with highest average distance
% prune it and recompute mapping precision
% prune the second highest average (from original tips)

    if nargin<5
        fe=0 ;
        if nargin<4
            ce=0 ;
        end
    end
    
    figure; treeplot(tree) ; [x,y] = treelayout(tree) ; text(x,y,num2str([(1:length(tree))']));
    
    % pairwise distance between tip weights
    tips = setdiff(1:numel(tree),unique(tree)) ;
    distNN=zeros(numel(tips));
    for w=1:numel(tips)
        for w2=w+1:numel(tips)
            distNN(w,w2) = ( DTWaverage(weight{tips(w)},weight{tips(w2)}) ) ;
            distNN(w2,w) = distNN(w,w2) ;
        end
    end
    % MDS plot
    [Y] = cmdscale(distNN);
    plot(Y(:,1),Y(:,2),'bo');
    for n=1:size(Y,1), text(Y(n,1)+.001,Y(n,2),num2str(tips(n))) ; end
    % density
    dens = zeros(1,numel(tips));
    cutoff = mean( reshape(distNN(distNN>0),1,numel(distNN(distNN>0))) ) ;
    for w=1:numel(tips)
        dens(w) = sum(distNN(w,distNN(w,:)<=cutoff)) ;
    end
    % minimum distance to a higher density weight
    d_dens = zeros(1,numel(tips));
    for w=1:numel(tips)
        min_dist_higher_dens = min(distNN(w,dens>dens(w))) ;
        if numel(min_dist_higher_dens)==0
            d_dens(w) = max(distNN(w,:)) ;
        else
            d_dens(w) = min(distNN(w,dens>dens(w))) ;
        end
    end
    %[tips; dens; d_dens; dens.*d_dens]
    % decision graph
    %plot(dens,d_dens,'o') ;
    % dens x d_dens in decreasing order plot
    %plot(sort(dens.*d_dens,'descend'),'o') ;
    centroids = tips( (dens.*d_dens)>mean(dens.*d_dens) ) ;
    % prune tree and centroids
    for c=centroids
        sibling = find(tree==tree(c)) ;
        sibling = sibling(sibling~=c) ;
        if sum(centroids==sibling) % sibling of centroid needs to be deleted, we need to copy centroid's weight to the parent and delete centroid node
            weight(tree(c)) = weight(c) ;
            tree(c) = NA ;
        end
    end
    tree(setdiff(centroids,tips)) = NA ; % delete all nodes that are not centroids
    change=1;
    while change==1
        change=0; 
        for t=tree
            if tree(t)~=NA && tree(tree(t))==NA
                tree(t)=NA ;
                change=1;
            end
        end
    end
    pruned_weight = weight(centroids) ;
    pruned_tree = tree ;
    
    

    %% using cluster_dp.m
    ND=8;
    N=sum(sum(distNN>0)); % number of pairwise comparisons
    dist=distNN;
    percent=2.0;
    position=round(N*percent/100);
    sda=sort(reshape(distNN(distNN>0),1,N));
    dc=sda(position);
    for i=1:ND
      rho(i)=0.;
    end
    %
    % Gaussian kernel
    %
    for i=1:ND-1
      for j=i+1:ND
         rho(i)=rho(i)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
         rho(j)=rho(j)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
      end
    end
    %
    % "Cut off" kernel
    %
    %for i=1:ND-1
    %  for j=i+1:ND
    %    if (dist(i,j)<dc)
    %       rho(i)=rho(i)+1.;
    %       rho(j)=rho(j)+1.;
    %    end
    %  end
    %end
    

    %% get BMU and distance for each read
    BMU = zeros(numel(reads),1,2) ; % full path to BMU + distance
    mp = zeros(numel(reads),2) ; % only the BMU + distance
    for s=1:numel(reads)
        tmp = findBMU(reads(s).seqvect, tree, weight, ce, fe, 'full') ;
        BMU(s,:,:) = 0 ;  % re-initialised BMU(s,:,:) if the new BMU in the tree is not as deep as the previous one
        BMU(s,1:size(tmp,1),1:2) = tmp ;
        mp(s,:) = BMU(s,size(tmp,1),1:2) ;
    end
    
    % calculate mapping precision
    map_prec_curr = sum(mp(:,2))/numel(reads) ;
    
    % order the tips
    % find mean distance for each weight (grouping reads according to their BMU)
    [mean_dist,weight_num] = grpstats(mp(:,2),mp(:,1),{'mean','gname'}) ;
    [~,I]=sort(mean_dist,'descend') ;
    bmu_ordered = str2double(weight_num(I)) ;
    
    for n=1:numel(bmu_ordered) 
        bmu_del=bmu_ordered(n) ;
        %fprintf(1,'%d\n',bmu_del) ;
       
        % find depth of the BMU tip in the tree
        bmu_del_depth = 1 ; % depth 1 is the root
        a = bmu_del ;
        if (tree(a)==0) % root is reached
            break;
        end
        while ( tree(a)~=0 ) 
            bmu_del_depth=bmu_del_depth+1; 
            a=tree(a); 
        end

        % prune the Etree, delete the node and its sister
        mapped_reads = find( ismember(mp(:,1), find(tree==tree(bmu_del))) ) ; % find the tree tips with same parent as maxid
        if numel(mapped_reads)==0
            continue; % no reads classified in this BMU because the sister has been deleted previously
        end
        BMU(mapped_reads,bmu_del_depth,1) = 0 ; % delete the tips: go up one level in the path
        BMU(mapped_reads,bmu_del_depth,2) = 0 ;

        % update mp for the reads that were mapping to the deleted BMU
        for m = 1:numel(mapped_reads)
            dep = BMU(mapped_reads(m),:,1)==max(BMU(mapped_reads(m),:,1)) ;
            mp(mapped_reads(m),:) = BMU(mapped_reads(m),dep,1:2) ;
        end
        map_prec_new = sum(mp(:,2))/numel(reads) ;
        
        fprintf(1,'%f -> %f (%f)\n',map_prec_curr,map_prec_new,map_prec_new-map_prec_curr) ;
        map_prec_curr = map_prec_new ;
    end
    
end
