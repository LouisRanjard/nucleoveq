function [ltree] = vqreads(truefasta,filefq,outputf)
% 
%% matlab -nojvm -r "addpath('/home/louis/Documents/Matlab/mfiles/dtw_nucleotide'); vqreads('haplo.fasta','file.fq','reconstructed.fasta'); quit();"
% 
% dirm='/home/louis/Documents/Matlab/mfiles/';
% cellfun(@(x) addpath(x),{dirm fullfile(dirm,'/DTWaverage') fullfile(dirm,'/dtw_nucleotide')});
% truefasta='/home/louis/Documents/Projects/ShortReads/test_popsize_nci/1PopDNA_10_1200_5k_0.haplo.fasta';
% filefq='/home/louis/Documents/Projects/ShortReads/test_popsize_nci/1PopDNA_10_1200_5k_0.combined.fq';
% 

true_seq = fastaread(truefasta);
for m=1:numel(true_seq) 
    true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; 
    true_seq(m).varcov = nan ;
end
reads = fastqread(filefq);
for m=1:numel(reads) 
    reads(m).seqvect=nucleo2mat(reads(m).Sequence) ; 
    reads(m).varcov = nan ;
end

% generate reference sequence for the root, starting by mutating one of the
reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
reference.seqvect = mutatematseq(true_seq(1).seqvect,0.3) ;
reference.Sequence = mat2nucleo(reference.seqvect) ;

% root weight initialisation: align all reads once to the reference to create initial weight
w = 0.0001^(1/numel(reads)) ;
rindex = randperm(numel(reads)) ;
reference.seqvect = [ reference.seqvect; ones(1,size(reference.seqvect,2)) ] ; % add persistence vector
for r=rindex % weight w must be enough to pull a position toward a different nucleotide
    % first generate a new read sequence that match the alignment length to the reference, by setting very low weight to the reference
    [ ~, new_read_seqvect ] = DTWaverage( reference.seqvect, reads(r).seqvect, 1, 0.001, 0, 1 ) ;
    % update the reference sequence with the read
    [ dist, reference.seqvect, aligned_pos ] = DTWaverage( reference.seqvect, reads(r).seqvect, 1, w, 0, 1 ) ;
    reads(r).bmu = 1 ;
    reads(r).distance = dist ;
    reads(r).position = aligned_pos ;
    %fprintf(1,'%f\n',DTWaverage(reads(r).seqvect,new_read_seqvect( 1:4, reads(r).position:min([size(new_read_seqvect,2) (reads(r).position+size(reads(r).seqvect,2)-1)])),1, 0.001, 0, 1) ) ;
    reads(r).seqvect = new_read_seqvect( 1:4, reads(r).position:min([size(new_read_seqvect,2) (reads(r).position+size(reads(r).seqvect,2)-1)]) ) ;
end
reference.Sequence = mat2nucleo(reference.seqvect) ;
[ ~, varcov ] = wcoverage(reference.seqvect(1:4,:),50,[reads.position],false) ;
reference.varcov = varcov ;

treer = lnode(nan,reference,1) ;
treer.accord() ;
ltree = treer ;
ltree = vqreadsplit(treer,reads,treer.entropy,ltree) ;

% inspect reconstructed sequences
weight2fasta([],arrayfun(@(x) {x.weight.seqvect},ltree),[outputf]) ;
seqalignviewer(multialign([ltree.weight])) ;
phytree = seqlinkage(seqpdist([[ltree.weight]'; true_seq]),'average',[[ltree.weight]'; true_seq]); phytreeviewer(phytree);

fprintf(1,'\n-- vqreads.m finishes --\n') ;

end
