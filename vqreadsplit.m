function ltree = vqreadsplit(parent,reads,cutoff_e,ltree)
%
% 
    % create two daughter nodes (if needed)
    fprintf(1,'Entropy: %f - %f\n',parent.entropy,cutoff_e);
    if parent.entropy<(cutoff_e/2)
        return ;
    else
        dleft = lnode(parent,parent.weight,length(ltree)+1) ;
        ltree = [ltree dleft] ;
        dright = lnode(parent,parent.weight,length(ltree)+1) ;
        ltree = [ltree dright] ;
        % Assign the reads to one or the other depending on minimum distance
        w = 0.0001^(1/numel(reads)) ;
        rindex = randperm(numel(reads)) ;
        readsleft=[];
        readsright=[];
        for r=rindex
            if vqdist(dleft,reads(r))<vqdist(dright,reads(r))
                dleft.update(reads(r),w) ;
                reads(r).bmu = dleft ;
                readsleft = [readsleft r];
            else
                dright.update(reads(r),w) ;
                reads(r).bmu = dright ;
                readsright = [readsright r];
            end
        end
        % Update fields for making them consistent to each other
        dleft.accord(reads(readsleft)) ;
        dright.accord(reads(readsright)) ;
        % Recursive call for each daughter node
        if (numel(readsleft)>0 && numel(readsright)>0)
            ltree = vqreadsplit(dleft,reads(readsleft),parent.entropy,ltree) ; %cell2mat( arrayfun(@(x) x.idl==3, [reads.bmu], 'Uniform', 0) ) ;
            ltree = vqreadsplit(dright,reads(readsright),parent.entropy,ltree) ;
        end
    end
    
end
