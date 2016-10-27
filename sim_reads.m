function [ reads ] = sim_reads(haplo,readlen,coverage)
% generate error-free short read sequences from a set of haplotypes with average
% coverage.
% header of read sequences names is "readX_haploY_posZ" where reaad number X from seuqence Y starting at position Z 
%
% haplo = haplotype sequences
% readlen = short read length
% ex: sim_reads(haplo,150,30)

    nhaplo = length(haplo) ; 
    len = abs(mean(arrayfun(@(x) length(x.Sequence), haplo))) ;
    
    readnum = (nhaplo * len * coverage) / readlen ; % 30x coverage
    for m=1:readnum
        haplo_id=randi([1 nhaplo],1,1) ;
        startpos=randi([1 length(haplo(haplo_id).Sequence)-readlen+1],1,1);
        endpos=startpos+readlen-1;
        reads(m).Header=['read' num2str(m) '_haplo' num2str(haplo_id) '_pos' num2str(startpos)] ;
        reads(m).seqvect=haplo(haplo_id).seqvect(:,startpos:endpos);
    end

end