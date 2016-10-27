function [ centroids ] = densipeak(dmat)
% cluster by density of peaks
% Rodriguez, A., & Laio, A. (2014). Clustering by fast search and find of density peaks. Science, 344(6191), 1492–1496.
%
% input: square distance matrix dmat

    % number of data points
    n = size(dmat,1) ;
    
    % local density
    dens = zeros(1,n);
    %cutoff = mean( reshape(dmat(dmat>0),1,numel(dmat(dmat>0))) ) ;
    cutoff = median( reshape(dmat(dmat>0),1,numel(dmat(dmat>0))) ) ;
    for w=1:n
        dens(w) = sum(dmat(w,dmat(w,:)<=cutoff)) ;
    end
    
    % minimum distance to a higher density weight
    d_dens = zeros(1,n);
    for w=1:n
        min_dist_higher_dens = min(dmat(w,dens>dens(w))) ;
        if numel(min_dist_higher_dens)==0
            d_dens(w) = max(dmat(w,:)) ; % convention for point with highest density
        else
            d_dens(w) = min(dmat(w,dens>dens(w))) ;
        end
    end
    %[1:n; dens; d_dens; dens.*d_dens]
    
    % decision graph
    %plot(dens,d_dens,'o') ;
    
    % dens x d_dens in decreasing order plot
    %plot(sort(dens.*d_dens,'descend'),'o') ;
    
    centroids = tips( (dens.*d_dens)>mean(dens.*d_dens) ) ;

end