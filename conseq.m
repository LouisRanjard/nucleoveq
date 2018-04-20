function [ reconsensus ] = conseq(SeqsMultiAligned, refindex)
% generate consensus of a multiple sequence alignment by taking the most common character at each position
% if redindex is given, then refindex if considered as a reference and a
% special case can occur: if only the reference is defined, ie all other
% sequences are '-' (ascii 45) then the reference nucleotide is retained in
% the consensus
%

try refindex=refindex; catch refindex=0 ; end

conse = SeqsMultiAligned(1).Sequence ;
for i = 1:size(SeqsMultiAligned(1).Sequence,2)
    nuc = cell(1,size(SeqsMultiAligned,1)) ;
    for j = 1:size(SeqsMultiAligned,1)
        nuc{j} = SeqsMultiAligned(j).Sequence(i) ;
    end
    v = +[nuc{:}]; % convert characters to ascii code
    if (refindex~=0 &&... % is only the reference defined at this position?
            v(refindex)~=45 &&...
            sum(v([1:refindex-1 refindex+1:end])==45)==size(SeqsMultiAligned,1)-1)
        conse(i) = char(v(refindex)) ;
    elseif (refindex~=0) % take the majority rule nucleotide of reads only
        conse(i) = char(mode(v([1:refindex-1 refindex+1:end]))) ;
    else % take the majority rule nucleotide
        conse(i) = char(mode(v)) ;
    end
end
%reconsensus.Sequence = erase(conse,'-') ; % >=R2016b
reconsensus.Sequence = conse(regexp(conse,'[ACTGactg]')) ;
reconsensus.Header = 'reconstructed consensus' ;
