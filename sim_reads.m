function [ reads, startpos ] = sim_reads(haplo,readlen,coverage,errorrate)
% generate error-free short read sequences from a set of haplotypes with average
% coverage.
% Header of read sequences names is "readX_haploY_posZ" where read number X from sequence Y starting at position Z 
% find original position with: str2double(reads(1).Header(strfind(reads(1).Header,'_pos')+4:end))
% errorrate: error rate in %, supposed uniform (Illumina: 0.1-0.5%)
%
% haplo = haplotype sequences
% readlen = short read length
% ex: sim_reads(haplo,150,30)

    if nargin<4, errorrate=0; end
    
    nhaplo = length(haplo) ; 
    len = abs(mean(arrayfun(@(x) length(x.Sequence), haplo))) ;
    
    readnum = (nhaplo * len * coverage) / readlen ; % coverage
    startpos = zeros(1,readnum) ;
    for m=1:readnum
        haplo_id=randi([1 nhaplo],1,1) ;
        startpos(m)=randi([1 length(haplo(haplo_id).Sequence)-readlen+1],1,1);
        endpos=startpos(m)+readlen-1;
        reads(m).Header=['read' num2str(m) '_haplo' num2str(haplo_id) '_pos' num2str(startpos(m))] ;
        if ( isfield(haplo(haplo_id),'seqvect') )
          reads(m).seqvect=haplo(haplo_id).seqvect(1:4,startpos(m):endpos);
        else
          reads(m).seqvect=nucleo2mat( haplo(haplo_id).Sequence(startpos(m):endpos) );  
        end
        if errorrate>0
            nerr=round(normrnd(readlen*errorrate/100,1)); % sample the number of error(s) in that read
            if nerr>0
                ids=datasample(1:readlen,nerr,'Replace',false) ; % find location of error (should be biased toward the end for Illumina but Uniform here)
                for n=ids
                    zero = find(reads(m).seqvect(1:4,n)==0) ; % row positions that ==0 and can be changed
                    reads(m).seqvect(1:4,n) = zeros(4,1) ; % set all positions to zero
                    k = randi([1 3],1,1) ;
                    reads(m).seqvect(zero(k),n) = 1 ;
                end
            end
        end
        reads(m).Sequence=mat2nucleo(reads(m).seqvect);
    end

end