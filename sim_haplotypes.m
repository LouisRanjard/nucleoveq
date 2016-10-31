function [ haplo, reference ] = sim_haplotypes(nhaplo,len,dist2ref)
% generate haplotypes sequences of specific distance to a random
% reference sequence
%
% nhaplo = number of haplotypes
% len = sequence genome length
% dist2ref = distance to the reference sequence (number of mutations/length)
% ex: sim_haplotypes(5,1200,0.1) for 90% similarity 

    reference = randmatseq(len) ;
    haplo = struct() ;
    for n=1:nhaplo % the order of the structure fields is important and must be consistent with other datasets (.Header,.Sequence,.seqvect,...)
        haplo(n).Header = ['haplotype' num2str(n)] ;
        vector_seq = mutatematseq(reference, dist2ref) ;
        haplo(n).Sequence = mat2nucleo(vector_seq) ;
        haplo(n).seqvect = vector_seq ;
    end
    
end
