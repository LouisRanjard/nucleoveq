function [tdist tdist2] = TreeDist(tree,a,b,distNN)
    % tdist: find the number of nodes between nodes a and b in the Tree tree
    % tdist2: sum the distances between these nodes
    
    if nargin<4, distNN=[]; end
    
    % intersection between the path back to the root for each node
    baca=a ;
    bacb=b ;
    while ( tree(a)~=0 ) baca=[baca tree(a)]; a=tree(a); end
    while ( tree(b)~=0 ) bacb=[bacb tree(b)]; b=tree(b); end
    inter = intersect(baca,bacb) ;
    lca = max(inter) ; % lowest common ancestor
    tdist = max( sum(baca>lca)+sum(bacb>lca)-1 , 0 ) ; % minimal distance must be 1, so 1 is subtracted (Samarasinghe2006)
    
    if numel(distNN)>0
        tdist2 = 0 ;
        for n=1:numel(baca)-1, tdist2 = tdist2+distNN(baca(n),baca(n+1)) ; end
        for n=1:numel(bacb)-1, tdist2 = tdist2+distNN(bacb(n),bacb(n+1)) ; end
    end

end
