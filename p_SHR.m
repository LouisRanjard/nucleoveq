function [ lproba_SHR ] = p_SHR(reads, tree, weight, ce, fe, verbose)
%
% 

    %% Proba of H (set of weight) given R (set of reads)
    BMU=zeros(numel(reads),2);
    for s=1:numel(reads)
        treepath = findBMU(reads(s).seqvect, tree, weight, ce, fe) ;
        BMU(s,:) = treepath ;
    end
    lproba_HgR = zeros(numel(unique(BMU(:,1))),1) ;
    n=1;
    for t=unique(BMU(:,1))'
        lproba_HgR(n) = sum(log(1-BMU(BMU(:,1)==t,2))) ; % product of 1 - all read alignment distances (1-d) (in log)
        n=n+1;
    end

    %% Proba of S (set of sequences) given H and R
    lproba_SgHR = zeros(numel(unique(BMU(:,1))),1) ;
    n=1;
    for t=unique(BMU(:,1))'
        lproba_SgHR(n) = sum(log(max(weight{t}(1:4,:),[],1))) ; % product of maximum probability of each position (in log)
    end

    %% Proba of S and H and R
    lproba_SHR = lproba_SgHR + lproba_HgR ;

end
