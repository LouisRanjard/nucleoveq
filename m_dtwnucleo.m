function [] = m_dtwnucleo(truefasta,filefq,outputf)
% 
% Ensemble version of dtwnucleo
%
% m_dtwnucleo('~/Documents/Projects/ShortReads/test_popsize/haplo.fasta','~/Documents/Projects/ShortReads/test_popsize/reads.fq','~/Documents/Projects/ShortReads/test_popsize/rechaplo.fasta')
% matlab -nojvm -r "addpath('.'); dtwnucleo('haplo.fasta','file.fq','reconstructed.fasta'); quit();"
%

%dirm='~/Documents/Matlab/mfiles' ;
%cellfun(@(x) addpath(x),{dirm fullfile(dirm,'/DTWaverage') fullfile(dirm,'/dtw_nucleotide')});
%truefasta='/home/louis/Documents/Projects/ShortReads/test_popsize_nci/test/2151_20000_5.haplo.fasta';
%filefq='/home/louis/Documents/Projects/ShortReads/test_popsize_nci/test/2151_20000_5.combined.fq';
true_seq = fastaread(truefasta);
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
reads=fastqread(filefq);
for m=1:numel(reads), reads(m).seqvect=nucleo2mat(reads(m).Sequence) ; end

% generate reference sequence for the root, starting by mutating one of the
reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
reference.seqvect = mutatematseq(true_seq(1).seqvect,0.3) ;
reference.Sequence = mat2nucleo(reference.seqvect) ;
ref_backup=reference;

% reads pairwise distance matrix defined as number of co-classification
read_dist = zeros(numel(reads)) ;

%% Generate multiple sets of candidate haplotypes
weight_seq = [] ;
s=1;
for g=1:10
    reference=ref_backup;
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
    [ ~, varcov ] = wcoverage(reference.seqvect(1:4,:),size(reads(1).seqvect,2),rpositionr,false) ;
    reference.varcov = varcov;
    reference.entropy = shannonEntropy(reference.seqvect) ;
    [tree, weight, BMU, ~, ~, ~] = ETDTWrec(reads,10,[1-w 1-w],0.9,[2 2],numel(reads),0.95,'seqvect',[],0,0,1,reference) ;
    
    tips = setdiff(1:length(tree),tree) ;
    for n=tips
        weight_seq(s).Header = ['weight_' num2str(g) '_' num2str(s)] ;
        weight_seq(s).seqvect = weight{n} ;
        weight_seq(s).Sequence = mat2nucleo(weight{n}) ;
        s=s+1;
    end
    
    % record reads co-classification (how many times each pair of reads are classified together)
    read_dist = read_dist + bsxfun(@eq,BMU(:,1)',BMU(:,1)) ;
    
end
% visualize weight sequences
phytree = seqlinkage(seqpdist([weight_seq'; true_seq; ref_backup]),'average',[weight_seq'; true_seq; ref_backup]); phytreeviewer(phytree);

%% Clean up, remove high entropy regions from the sequence matrices
for w=1:numel(weight_seq)
    [ meanEntropy, shannonEnt ] = shannonEntropy(weight_seq(w).seqvect(1:4,:)) ;
    figure; plot(shannonEnt); title(meanEntropy);
end

%% Cluster the candidate haplotypes
reference = ref_backup ;
[ ~, varcov ] = wcoverage(reference.seqvect(1:4,:),50,ones(1,length(weight_seq)),false) ;
reference.varcov = varcov ;
w = 0.0001^(1/numel(weight_seq)) ;
[htree, hweight, ~, ~, ~, ~] = ETDTWrec(weight_seq,5,[1-w 1-w],0.9,[2 2],numel(weight_seq),0.95,'seqvect',[],0,0,1,reference) ;
haplo_seq = weight2fasta(htree,hweight) ;
phytree = seqlinkage(seqpdist([haplo_seq'; true_seq; ref_backup]),'average',[haplo_seq'; true_seq; ref_backup]); phytreeviewer(phytree);

fprintf(1,'\n-- m_dtwnucleo.m finishes --\n') ;

end
