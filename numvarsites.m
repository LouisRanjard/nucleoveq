function [vs] = numvarsites(haplo)
% return the variable sites positions for a set of sequences (e.g.
% haplotypes) of identical length

    hlen = length(haplo(1).Sequence) ;
    nhaplo = numel(haplo) ;
    
    vs=zeros(1,hlen);
    for n=1:hlen
        tmp = zeros(4,1) ;
        for h=1:nhaplo
            tmp = tmp + haplo(h).seqvect(1:4,n) ;
        end
        if sum(tmp(tmp>0)~=h)
            vs(n) = n ;
        end
    end
    
    fprintf(1,'Number of variable sites: %d/%d (%.2f%%)\n',length(vs(vs>0)),hlen,100*length(vs(vs>0))/hlen);
    
end
