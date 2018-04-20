function [ ] = plotnucleoveq(mseq,reference)
% plot coverage, persistence and entropy of a matrix nucleotide sequence

[ ~, shannonEnt ] = shannonEntropy(mseq.seqvect);

%figure ;
clf('reset');
subplot(2,1,1);
plot(mseq.seqvect(5,:)) ; %ylim([0.98 1.02]) ;
ylabel('Persistence');
yyaxis right ; 
plot(mseq.seqvect(6,:)) ;
ylabel('Coverage');

subplot(2,1,2);
plot(shannonEnt); 
ylabel('Entropy');

if nargin>1 % reference sequence is input, plot alignment using bioinformatics toolbox
    [~,align] = nwalign(reference, rmfield(mseq,{'coverage','seqvect'}),'Glocal',true) ;
    seqalignviewer(align,'SeqHeaders',{reference.Header,mseq.Header});
end
