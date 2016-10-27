function [] = dtwnucleo(truefasta,filefq,outputf)
% 
% dtwnucleo('~/Documents/Projects/ShortReads/test_popsize/haplo.fasta','~/Documents/Projects/ShortReads/test_popsize/reads.fq','~/Documents/Projects/ShortReads/test_popsize/rechaplo.fasta')
% matlab -nojvm -r "addpath('.'); dtwnucleo('haplo.fasta','file.fq','reconstructed.fasta'); quit();"
%

%cellfun(@(x) addpath(x),{dirm fullfile(dirm,'/DTWaverage') fullfile(dirm,'/dtw_nucleotide')});
true_seq = fastaread(truefasta);
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
reads=fastqread(filefq);
for m=1:numel(reads), reads(m).seqvect=nucleo2mat(reads(m).Sequence) ; end

% generate reference sequence for the root, starting by mutating one of the
reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
reference.seqvect = mutatematseq(true_seq(1).seqvect,0.3) ;
reference.Sequence = mat2nucleo(reference.seqvect) ;

% root weight initialisation: align all reads once to the reference to create initial weight
w = 0.0001^(1/numel(reads)) ;
rindex = randperm(numel(reads)) ;
reference.seqvect = [ reference.seqvect; ones(1,size(reference.seqvect,2)) ] ; % add persistence vector
BMUr = zeros(numel(reads),2) ;
rpositionr = zeros(numel(reads),1) ;
for r=rindex % weight w must be enough to pull a position toward a different nucleotide 
    [ dist, reference.seqvect, aligned_pos ] = DTWaverage( reference.seqvect, reads(r).seqvect, 1, w, 0, 1 ) ;
    rpositionr(r) = aligned_pos ;
    BMUr(r,:) = [1 dist] ;
end
reference.Sequence = mat2nucleo(reference.seqvect) ;
[ ~, varcov ] = wcoverage(reference.seqvect(1:4,:),50,rpositionr,false) ;
reference.varcov = varcov;
reference.entropy = shannonEntropy(reference.seqvect) ;
[tree, weight, ~, ~, ~, ~] = ETDTWrec(reads,50,[1-w 1-w],0.9,[2 2],numel(reads),0.95,'seqvect',[],0,0,1,reference) ;

weight2fasta(tree,weight,[outputf]) ;

fprintf(1,'\n-- dtwnucleo.m finishes --\n') ;

end

