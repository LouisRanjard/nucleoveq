function [ weight, BMU, readposition ] = fortify(ite, reads, weight, w, BMU, ce, fe)
% realign each read to its BMU "ite" iteration(s)
% stop if ....

    if nargin<5
        fe=0 ;
        if nargin<4
            ce=0 ;
        end
    end
    
    % record alignment position of the reads
    rposition = zeros(numel(reads),ite) ;
    
    % align the reads to their respective BMU
    for n=1:ite
        sindex = randperm(numel(reads)) ;
        for s=sindex
            [ BMU(s,2), weight{BMU(s,1)}, aligned_pos ] = DTWaverage( weight{BMU(s,1)}, reads(s).seqvect, 1, w, ce, fe ) ;
            rposition(s,n) = aligned_pos ;
        end
    end
    
    % get a single alignment position for each read as the most frequent position
    C = num2cell(rposition,2) ;
    readposition = cellfun(@mode,C);
    
end