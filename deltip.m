function [tree,weight,locked] = deltip(tree,weight,dtip,locked)
% delete the tip dtip and copy the weight matrix of its sibling to replace parent node weight matrix
% locked is a set of nodes that cannot be deleted
    
    % control
    if sum(locked==dtip)>0, error('Cannot delete locked node %d.',dtip); end
    if sum(tree==dtip)>0, error('Node %d is not a tip.',dtip); end
    
    % find sibling
    dsibling = find(tree==tree(dtip)) ;
    dsibling = dsibling(dsibling~=dtip) ;
    
    % delete the node
    if dtip<numel(tree)
        tree(tree>dtip)=tree(tree>dtip)-1; % update index of node that are changing indexes
        for ti=dtip:(numel(tree)-1)
            tree(ti)=tree(ti+1);
            weight(ti)=weight(ti+1);
        end
        locked(locked>dtip)=locked(locked>dtip)-1;
    end
    tree = tree(1:(end-1)) ;
    weight = weight(1:(end-1)) ;
    
    % copy sibling weight matrix to parent
    if numel(dsibling)>0
        if dsibling>dtip % index has changed since dtip was deleted
            dsibling=dsibling-1;
        end
        if sum(tree==dsibling)==0 % check if sibling has children
            if sum(locked==dsibling)==0
                weight(tree(dtip-1)) = weight(dsibling) ;
                [tree,weight,locked] = deltip(tree,weight,dsibling,locked) ;
            end
        end
%         if sum(tree==dsibling)>0 % check if sibling has children
%             %for da=find(tree==dsibling)
%             %    if sum(locked==da)==0 && sum(tree==da)==0
%             %        [tree,weight] = deltip(tree,weight,da,locked) ;
%             %    end
%             %end
%         else % copy sibling weight matrix to parent
%             if numel(dsibling)>0
%                 weight(tree(dtip-1)) = weight(dsibling) ;
%             end
%         end
    end
   
    % delete nodes with only one child
    for u=unique(tree)
        if u~=0 && sum(tree==u)==1
            child = find(tree==u) ;
            children = find(tree==child) ; % find children of child
            tree(child)=tree(u) ; % replace node by father
            weight(u) = weight(child) ;
            tree(tree>u) = tree(tree>u)-1 ;
            tree(children) = u ; % update children with new father index
            for ti=child:(numel(tree)-1)
                tree(ti)=tree(ti+1);
                weight(ti)=weight(ti+1);
            end
            tree = tree(1:(end-1)) ;
            weight = weight(1:(end-1)) ;
            locked(locked>u)=locked(locked>u)-1 ;
        end
    end
    
end
