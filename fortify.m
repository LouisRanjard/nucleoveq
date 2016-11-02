function [ weight, BMU, readposition ] = fortify(ite, limit, reads, weight, w, BMU, ce, fe, verbose)
% realign each read to its BMU "ite" iteration(s)
% stop if ....

    if ite==0 % default maximum number of iterations
        ite=1000 ;
    end
    
    if limit==0 % default entropy limit is set to maximum possible for 4 rows : -sum([.25 .25 .25 .25].*log([.25 .25 .25 .25]))=1.3863
        limit=1.3863 ;
    end

    % verbose mode, plot entropy values
    if nargin<9
        verbose=0 ;
        if nargin<8
            fe=0 ;
            if nargin<7
                ce=0 ;
            end
        end
    end
    
    % record alignment position of the reads
    rposition = zeros(numel(reads),ite) ;
    
    % counter, if the entropy does not change for three iterations, then exit
    cnt3 = 0 ;
    
    % align the reads to their respective BMU
    x=zeros(1,ite);
    if verbose
        fprintf(1,'Consensus entropy (iteration=%d, limit=%f)\n',ite,limit);
    end
    for n=1:ite
        sindex = randperm(numel(reads)) ;
        for s=sindex
            [ BMU(s,2), weight{BMU(s,1)}, aligned_pos ] = DTWaverage( weight{BMU(s,1)}, reads(s).seqvect, 1, w, ce, fe ) ;
            rposition(s,n) = aligned_pos ;
        end
        % parallelising the previosu loop by grouping reads by BMU
%         for b = unique(BMU(:,1))
%            for y = (BMU(s,1)==b)
%               [ BMU(y,2), weight{BMU(y,1)}, aligned_pos ] = DTWaverage( weight{BMU(y,1)}, reads(y).seqvect, 1, w, ce, fe ) ;
%               rposition(y,n) = aligned_pos ;
%            end
%         end
        ses = shannonEntropy_s(weight(unique(BMU(:,1)))) ; % only consider tips (weight that are BMU to some read)
        if (n>1 && ses==x(n-1))
            cnt3 = cnt3+1 ;
        else
            cnt3 = 0 ;
        end
        if verbose
            fprintf(1,'%d - %f\n',n,ses) ;
        end
        x(n)=ses;
        if ses>limit
            break;
        end
        if cnt3>2
            break;
        end
    end
    
    % get a single alignment position for each read as the most frequent position
    C = num2cell(rposition,2) ;
    readposition = cellfun(@mode,C);
    
    if verbose
        figure;
        plot(x(1:n));
        line([1 n],[limit limit],'Color','red');
    end
    
end