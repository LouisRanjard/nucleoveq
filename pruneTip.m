function [tree,weight,nbmu,locked] = pruneTip(tree,weight,nbmu,locked)
% remove all the tips in the tree that are not in [locked]
% then remove single child internal nodes by copying node below
%

    % control
    if numel(locked)<2, error('Only one tip is to be kept (need at least two).'); end
    
    % delete the tips
    tips = setdiff(1:numel(tree),unique(tree)) ;
    tobedeleted = setdiff(tips,locked) ;
    while numel(tobedeleted)>0
        [tree,weight,nbmu,locked] = delTip(tree,weight,nbmu,tobedeleted(1),locked) ;
        tips = setdiff(1:numel(tree),unique(tree)) ;
        tobedeleted = setdiff(tips,locked) ;
    end
    
    %figure; treeplot(tree); [x,y]=treelayout(tree);text(x,y,num2str([(1:length(tree))']));text(x+.05,y,num2str(cell2mat(weight(:))));
    
    % delete the nodes with only 1 child left
    node = unique(tree) ;
    nodecount = arrayfun(@(x) sum(tree==x),node) ;
    tobedeleted = node(nodecount==1) ; % always include the root as unique
    while numel(tobedeleted)>1
        [tree,weight,nbmu,locked] = delOnlyChild(tree,weight,nbmu,max(tobedeleted),locked) ; % delete the deepest node first
        node = unique(tree) ;
        nodecount = arrayfun(@(x) sum(tree==x),node) ;
        tobedeleted = node(nodecount==1) ;
    end

    %% internal functions
    
    function [tree,weight,nbmu,locked] = delTip(tree,weight,nbmu,dtip,locked)
        if dtip<numel(tree)
            tree(tree>dtip) = tree(tree>dtip)-1 ; % update index of node that are changing indexes
            for ti=dtip:(numel(tree)-1)
                tree(ti) = tree(ti+1) ;
                weight(ti) = weight(ti+1) ;
                nbmu(ti) = nbmu(ti+1) ;
            end
            locked(locked>dtip) = locked(locked>dtip)-1 ;
        end
        tree = tree(1:(end-1)) ;
        weight = weight(1:(end-1)) ;
        nbmu = nbmu(1:(end-1)) ;
    end

    function [tree,weight,nbmu,locked] = delOnlyChild(tree,weight,nbmu,dnode,locked)
        child = find(tree==dnode) ;
        children = find(tree==child) ; % find children of child
        tree(children) = dnode ; % update children with new father index
        if dnode>1 % exception: do not update root weight
            weight(dnode) = weight(child) ;
            nbmu(dnode) = nbmu(child) ;
        end
        tree(tree>child) = tree(tree>child)-1 ;
        for ti=child:(numel(tree)-1)
            tree(ti) = tree(ti+1);
            weight(ti) = weight(ti+1);
            nbmu(ti) = nbmu(ti+1);
        end
        tree = tree(1:(end-1)) ;
        weight = weight(1:(end-1)) ;
        nbmu = nbmu(1:(end-1)) ;
        locked(locked>child) = locked(locked>child)-1 ;
    end

end

%% tests tree pruning
% figure;
% t0=[0 1 1 2 2 3 3 4 4 6 6 8 8 12 12];
% w0={[1] [2] [3] [4] [5] [6] [7] [8] [9] [10] [11] [12] [13] [14] [15]};
% treeplot(t0);
% [x,y]=treelayout(t0);
% text(x,y,num2str([(1:length(t0))']));
% [t,w,l]=deltip(t0,w0,11,[10 7]);
% figure; 
% treeplot(t); 
% [x,y]=treelayout(t);
% text(x,y,num2str([(1:length(t))']));
% text(x+.05,y,num2str(cell2mat(w(:))));

% figure; treeplot(tree); [x,y]=treelayout(tree);text(x,y,num2str([(1:length(tree))']));text(x+.05,y,num2str(cell2mat(weight(:))));
% [t,w,l]=pruneTip([0 1 1 2 2 4 4],{[1] [2] [3] [4] [5] [6] [7]},[3 5 6]); figure; treeplot(t); [x,y]=treelayout(t);text(x,y,num2str([(1:length(t))']));text(x+.05,y,num2str(cell2mat(w(:))));
% [t,w,l]=pruneTip([0 1 1 2 2 4 4],{[1] [2] [3] [4] [5] [6] [7]},[3 5]); figure; treeplot(t); [x,y]=treelayout(t);text(x,y,num2str([(1:length(t))']));text(x+.05,y,num2str(cell2mat(w(:))));
% [t,w,l]=pruneTip([0 1 1 2 2 4 4],{[1] [2] [3] [4] [5] [6] [7]},[3 6 7]); figure; treeplot(t); [x,y]=treelayout(t);text(x,y,num2str([(1:length(t))']));text(x+.05,y,num2str(cell2mat(w(:))));
% [t,w,l]=pruneTip([0 1 1 2 2 3 3],{[1] [2] [3] [4] [5] [6] [7]},[5 6]); figure; treeplot(t); [x,y]=treelayout(t);text(x,y,num2str([(1:length(t))']));text(x+.05,y,num2str(cell2mat(w(:))));
% [t,w,l]=pruneTip([0 1 1 2 2 4 4],{[1] [2] [3] [4] [5] [6] [7]},[3 5 6]); figure; treeplot(t); [x,y]=treelayout(t);text(x,y,num2str([(1:length(t))']));text(x+.05,y,num2str(cell2mat(w(:))));
% [t,w,l]=pruneTip([0 1 1 2 2 4 4],{[1] [2] [3] [4] [5] [6] [7]},[3 5]); figure; treeplot(t); [x,y]=treelayout(t);text(x,y,num2str([(1:length(t))']));text(x+.05,y,num2str(cell2mat(w(:))));
% [t,w,l]=pruneTip([0 1 1 2 2 4 4],{[1] [2] [3] [4] [5] [6] [7]},[3 6 7]); figure; treeplot(t); [x,y]=treelayout(t);text(x,y,num2str([(1:length(t))']));text(x+.05,y,num2str(cell2mat(w(:))));
% [t,w,l]=pruneTip([0 1 1 2 2 4 4],{[1] [2] [3] [4] [5] [6] [7]},[3 5]); figure; treeplot(t); [x,y]=treelayout(t);text(x,y,num2str([(1:length(t))']));text(x+.05,y,num2str(cell2mat(w(:))));
% [t,w,l]=pruneTip([0 1 1 2 2 4 4],{[1] [2] [3] [4] [5] [6] [7]},[3 5 6]); figure; treeplot(t); [x,y]=treelayout(t);text(x,y,num2str([(1:length(t))']));text(x+.05,y,num2str(cell2mat(w(:))));
% [t,w,l]=pruneTip([0 1 1 2 2 4 4 6 6],{[1] [2] [3] [4] [5] [6] [7] [8] [9]},[3 5 8]); figure; treeplot(t); [x,y]=treelayout(t);text(x,y,num2str([(1:length(t))']));text(x+.05,y,num2str(cell2mat(w(:))));
% [t,w,l]=pruneTip([0 1 1 2 2 4 4 6 6],{[1] [2] [3] [4] [5] [6] [7] [8] [9]},[9 5 8]); figure; treeplot(t); [x,y]=treelayout(t);text(x,y,num2str([(1:length(t))']));text(x+.05,y,num2str(cell2mat(w(:))));
% [t,w,l]=pruneTip([0 1 1 2 2 4 4 6 6],{[1] [2] [3] [4] [5] [6] [7] [8] [9]},[9 5 ]); figure; treeplot(t); [x,y]=treelayout(t);text(x,y,num2str([(1:length(t))']));text(x+.05,y,num2str(cell2mat(w(:))));
% [t,w,l]=pruneTip([0 1 1 2 2 4 4 6 6],{[1] [2] [3] [4] [5] [6] [7] [8] [9]},[9 5 7 8]); figure; treeplot(t); [x,y]=treelayout(t);text(x,y,num2str([(1:length(t))']));text(x+.05,y,num2str(cell2mat(w(:))));
