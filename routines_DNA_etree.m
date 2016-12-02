% 
direc = dir('/home/louis/Sounds/Tieke/K_01_2007/part*') ;
for n=1:numel(direc)
    direcres = dir(['/home/louis/Sounds/Tieke/K_01_2007/' direc(n).Header '/result*']) ;
    tmp = ['/home/louis/Sounds/Tieke/K_01_2007/' direc(n).Header '/' direcres.Header ] ;
    fprintf(1,'%s\n',tmp) ;
    load(fullfile(tmp, 'dictionary.mat')) ;
    syllab = cell2struct( cellfun(@(x) norma(x,0,1),{syllab.seqvect}, 'UniformOutput', false) , 'seqvect' , 1 )' ;
    save( fullfile(tmp, 'syll.mat'), 'syllab') ;
    clear syllab dico ;
end

% load everything and save in a new structure
% WATCH OUT THE INDEX IN song
direc = dir('/home/louis/Sounds/Tieke/K_01_2007/part*') ;
sylb = [] ;
for n=1:numel(direc)
    direcres = dir(['/home/louis/Sounds/Tieke/K_01_2007/' direc(n).Header '/result*']) ;
    tmp = ['/home/louis/Sounds/Tieke/K_01_2007/' direc(n).Header '/' direcres.Header ] ;
    fprintf(1,'%s\n',tmp) ;
    load(fullfile(tmp, 'syll.mat')) ;
    sylb = [sylb syllab]; 
end

%% test with DNA sequences (Anolis dataset)

dnaseq=[];
dnaseq(1).Header='Anolis_ahli'  ;         dnaseq(1).Sequence= 'ATGAGCCCAATAATATACACAATTATACTATCAAGCCTAGCAACAGGCACTATCGTTACCATAACGAGCTACCACTGACTCCTAGCCTGAATCGGACTAGAAATAAACACTTTATCAATTATTCCAATTATTTCTACCATACACCACCCACGATCAACAGAGGCCGCCACAAAATACTTTTTAACACAAGCAGCGGCTTCTGCCATAATCTTGTTTTCAAGCATAATTAACGCCTGACAAACCGGATCAT';
dnaseq(2).Header='Anolis_aliniger'   ;    dnaseq(2).Sequence= 'ATGAGCCCTACAGTTTATTCAATTATTTTGTCAAGCCTACCAACAGGCACAGTTATTACTATAACCAGCTACCATTGATTAATAGCCTGAGTCGGGCTAGAAATTAACACACTCGCAATTATTCCTGTTGTTTCAATACAACATCACCCACGGTCCACAGAAGCCGCCACAAAATATTTTCTAACACAAGCAGCAGCCTCCGCCTTAATTCTATTTGCTAGCACAACTAACGCCTGATCAACGGGCACAT';
dnaseq(3).Header='Anolis_alutaceus'  ;    dnaseq(3).Sequence= 'ATGAACCCAACAATTATTATAATTACCCTAACCAGCCTGGCAACTGGTACAGTTATTACCATACATAGCTTCCATTGATTAATGGCCTGANTCGGATTAGANATCAATACACTATCAATTATTCCAATAATTTCAACATTACACCACCCACGATCAACTGAAGCTGCTACAAAATATTTCCTCACCCAAGCAGCTGCTTCANCTTTAATCCTTTTTTCAAGCACAATTAATGCCTGACAAACAGGATCAT';
dnaseq(4).Header='Anolis_angusticeps' ;   dnaseq(4).Sequence= 'ATGAGCCCCCCCATTTTTACAATTATCATCTCAAGTCTAGCAACAGGTACAATTATTACCATAACCAGCTACCATTGACTCATAGCCTGAGTTGGTCTAGAAATAAATACACTAGCAATTATTCCTATTATTTCAACAACACATCACCCACGAGCCACAGAAGCTTCCACAAAATATTTTCTTACACAAGCTGCAGCCTCTGCTCTAATTTTATTTTCTAGTATAATTAACGCATGACACACAGGATCTT';
dnaseq(5).Header='Anolis_bahorucoensis';  dnaseq(5).Sequence= 'ATGAGCCCCATAATTTACTCAATTGTATTCTCAAGCCTAGCNACAGGTACTATTATTACTATAACCAGCTACCACTGATTTATGGCCTGAATCGGACTAGAAATTAATACACTAGCAGTAATCCCCATTATTTCAACACTACACCACCCACGATCTACAGAAGCTGCTACAAAATACTTCTTAACACAAGCAGCAGCCTCCGCCACAATCCTATTTTCAAGTATAATTAATGCCTGACAAACAGGCACAT';
dnaseq(6).Header='Anolis_barahonae'   ;   dnaseq(6).Sequence= 'ATGAGCCCGCTAATTTATATGATTATTTTATCAAGCTTAGCAACAGGCACAATTATTACAATAACGAGTTTTCATTGAATTATAGCTTGAATTGGGTTAGAAATCAACACCTTAGCAATTATCCCAATTATTTCTATATTACACCACCCACGTTCTACTGAAGCAGCCACAAAATATTTTCTTACACAAGCAGCAGCATCCGCTATAATCCTATTTTCAAGTATAATTAATGCCTGACAAACAGGAACAT';
dnaseq(7).Header='Anolis_brevirostris';   dnaseq(7).Sequence= 'ATGAGCCCACTAATCCACACAATTATACTCTCAAGTCTAGCAACAGGCACTATTATTACTATATCTAGCCACCACTGACTAATAGCCTGAATTGGATTAGAAATTAACACACTAGCAATTATCCCCATCATTTCAACATCCCACCACCCACGATCAACAGAAGCTGCCACAAAATATTTCCTTACACAAGCAGCAGCCTCTGCCACCGTACTATTTTCTAGTATAATTAATGCCTGACAAACCGGAACAT';
dnaseq(8).Header='Anolis_coelestinus' ;   dnaseq(8).Sequence= 'ATGAGCCCACTAATTTTTTCAATCGTCCTGTCAAGCCTAGCAACAGGCACTATTATTACCATAACCAGCTATCACTGATTAATAGCTTGAATTGGTCTAGAAATAAACACACTTGCTATTATTCCAATTATCTCAATACAACATCACCCTCGATCTACAGAAGCCGCTACAAAATATTTCCTTACACAAGCAGCAGCCTCCGCTATGATTTTATTCGCCAGCACAACAAATGCTTGATACACAGGCACAT';
dnaseq(9).Header='Anolis_cristatellus' ;  dnaseq(9).Sequence= 'ATGAGCNNNACAATCTACACAATTATTTTGTNNNNNCTAGCAACAGGCACTATCATCACTATAACTAGCTTCCACTGACTAATGGCCTGAATCGGACTAGAGCTTAATACGCTAGCAATTATCCCGATTATTTCAACATTACACCACCCACGATCAACAGAAGCCGCAACAAAATACTTCTTAACACAAGCAGCAGCCTCTGCAATAATTATGTTTTCTAGCATAATTAATGCCTGAAACATAGGAACAT';
dnaseq(10).Header='Anolis_cuvieri'   ;    dnaseq(10).Sequence='ATGAGCCCAACAATTCTCTCAATCATTTTATCAAGCCTAGCAGCAGGAACAATTATTACAATAACAAGCTTTCATTGATTAATAGCCTGAATTGGACTAGAAATTAATACACTAGCAATTATTCCAATTATCTCAATAATACATCACCCACGATCTACAGAAGCAGCCACAAAATATTTTCTCACACAAGCAGCAGCATCAGCTATAATCCTGTTCTCAAGCATAATTAATGCTTGACAAACAGGGACAT';
dnaseq(11).Header='Anolis_distichus'  ;   dnaseq(11).Sequence='ATGAGCCCGCCAATCTACGCAATTATACTATCAAGCTTAGCAACAGGCACCATTATCACTATAACCAGTTACCATTGACTAATGGCCTGAATTGGACTAGAAATTAATACACTAGCAATTCTTCCAATTATTGCAACATCACACAACCCACGATCCACAGAAGCTGCCACAAAATACTTTTTAACACAATCAGCAGCTTCCGCCACTATCTTATTTTCTAGCATACTTAACGCCTGACAAACCGGAACAT';
dnaseq(12).Header='Anolis_equestris'  ;   dnaseq(12).Sequence='ATGAGCCCAACAATTTATTCAATTATCCTATCAAGCCTTGCAANNNNNACAATTATTACTATAACCAGCCACCATTGACTAATAGCCTGANTCGGATTAGAAATTAACACATTAGCAATTATCCCAATTATTTCAACATTACACCACCCACGATCCACAGAAGCCGCCACAAAATATTTCCTAACACAAGCAGCTGCTTCTGCTATAATTTTATTTTCTAGCATAACAAATGCTTGATACACAGGTACAT';
dnaseq(13).Header='Anolis_garmani'   ;    dnaseq(13).Sequence='ATGAGCCCAACCATCCTTATAATTATTATCTCAAGCCTGGCAACAGGTACCATTATTACCATAACAAGCCACCACTGACTCATAGCCTGAATCGGACTAGAAATAAATACCTTAGCTATTATTCCAATCATTACTACTATACACAACCCACGATCAACAGAAGCCGCCACAAAATACTTCTTAACACAAGCAGCAGCCTCTGCCATAATCTTATTCTCAAGCATAATTAATGCCTGACAAATAGGATCAT';
dnaseq(14).Header='Anolis_grahami'    ;   dnaseq(14).Sequence='ATGAGCCCATCAATCCTTATAATTATTATTTCAAGCCTGGCAACAGGCACTATTATTACTATAACAAGCCACCACTGACTTATAGCCTGAGTCGGACTAGAAATAAATACTTTGGCAATTATCCCAATTATTTCTACTACACACAGCCCGCGATCCACAGAAGCCGCTACAAAATATTTTTTAACACAAGCAGCTGCCTCTACCATAATCTTATTTTCAAGCATAACCAACGCCTGACAAACAGGCACAT';
dnaseq(15).Header='Anolis_insolitus'  ;   dnaseq(15).Sequence='ATGAACCCAACTATTCTTACATTAATTTTATCAAGCTTAGCAACAGGTACAATCCTTACAATAATCAGCTTTCACTGACTACTCGCATGAATTGGTCTAGAGATTAATACCCTAGCAATTATTCCTATTATCTCAGCACCTCACCACCCCCGACCCACAGAAGCCTCCACAAAATATTTCCTTACACAAGCAGCTGCTTCTGCTACAATCTTATTTTCAAGCATAATTAATGCTTGACTAACAGGCACAT';
dnaseq(16).Header='Anolis_krugi'      ;   dnaseq(16).Sequence='ATGAGCCCTGCAATTTACACAATTATATTATCAAGCTTGGCAACAGGCACTATCATCACTATAACAAGCTTCCACTGACTAATAGCCTGAGTTGGACTAGAACTTAATACATTAGCAATTATCCCAATTATTTCAACATTACACCACCCCCGGGCCACAGAAGCCTCAACAAAGTATTTCCTCACTCAAGCAGCAGCCTCTGCCATAATTTTATTTTCTAGCATAATTAATGCCTGACACACGGGAACAT';
dnaseq(17).Header='Anolis_lineatopus'  ;  dnaseq(17).Sequence='ATGAGCCCTGCTATCTTTTCAATTGTCATATCCAGCCTAGCAACAGGCACTATTATTACCATAACAAGCCACCACTGACTATTAGCCTGAATAGGGTTAGAAATAAATACCCTAGCAATTATTCCAATTATTTCTATTACCCACAACCCCCGAGCCACAGAAGCCGCCACAAAGTACTTCTTAACACAAGCAGCAGCCTCAGCCATAGTTTTATTTGCAAGCATAACTAATGCTTGACAAACAGGAACAT';
dnaseq(18).Header='Anolis_loysiana'   ;   dnaseq(18).Sequence='ATGAACCCAGTCGTTATTNNNATTCTTCTATCAAGCCTAGCAACTGGTACTATTATTACTATAACTAGCTATCACTGATTAATAGCCTGAGTTGGACTAGAAATAAACACATTAGCAATTATTCCAATCATTTCAACAACTCACCACCCACGAGCCACAGAAGCAGCCACAAAATACTTCTTAACCCAAGCTGCAGCCTCCGCCCTAATCTTATTTTCAAGTACAATTAATGCTTGATACTCAGGCTCAT';
dnaseq(19).Header='Anolis_luteogularis' ; dnaseq(19).Sequence='ATGAGTCCAACAATTTATTCAATTATTCTATCAAGCCTTGCAACAGGTACAATTATTACTATAACCAGCTACCATTGACTAATAGCCTGAGTCGGATTAGAAATTAACACATTAGCAATTATCCCAATTATTTCAACATTACACCATCCACGATCCACAGAAGCAGCCACAAAATACTTCCTAACACAAGCAGCTGCTTCTGCTATAATTTTATTTTCTAGCATAACAAATGCTTGATACACAGGTACAT';
dnaseq(20).Header='Anolis_marcanoi'   ;   dnaseq(20).Sequence='ATGAGCCCAACAATTTTTTCAATTATGCTATCGAGTCTAGCAACAGGCACCATTATTACTATAACAAGCTTTCACTGACTAATGGCCTGAGTCGGCTTAGAAATTAATACGCTAGCTGCTATTCCAATTATTTCAATACAACACCACCCTCGATCAACAGAAGCAGCCACAAAATACTTTTTAACACAAGCAACTGCATCCTCCTTAATTTTATTTTCAAGTATGATTAATGCTTGGCATACAGGAACAT';
dnaseq(21).Header='Anolis_occultus'   ;   dnaseq(21).Sequence='ATGAGCCCCAATAATCTACTTAATAGTTTAATTAGCTTATTTATANNNACAACACTAGTAACCACTAGCCACCACTGATTATTAGCGTGAGTTGGCTTGGAAATTAACACACTTGCAGCTATTCCACTTATCTCAACAAAACATCACCCCCGAGCTACAGAATCAGCCACAAAATACTTTTTAATTCAAGCAGCAGCCTCAGCTACAATCTTATTTTCAAGTACCATTAATGCTTGACACACAGGCTCAT';
dnaseq(22).Header='Anolis_olssoni'    ;   dnaseq(22).Sequence='ATGAACCCCACCATCTCCATAAATTATCTATCAAGCCTAGCAACAGGAACAATTATTACTATGACCAGCTTTCATTGATTAATAGCATGAATTGGATTAGAAGTCAACACACTAGCAATTATTCCAATCATCTCAGCCCCTCACCACCCACGATCAACAGAAGCTGCAACAAAATACTTTCTCACACAAGCAGCTGCCTCCGCTATAATTCTATTTGCCAGTATAATTAACGCCTGACAAACAGGCACAT';
dnaseq(23).Header='Anolis_ophiolepis'  ;  dnaseq(23).Sequence='ATGAGCCCAACAATCTTTATAATTATTTTATCAAGTCTTGCAACTGGTACAATTATTACTATAACTAGTTATCACTGACTATTAGCCTGAATCGGCCTAGAAATTAATACCTTATCAATTATCCCACTTATTTCAACAACCCACCATCCACGAGCCACAGAAGCCGCTACCAAGTATTTTCTTACACAAGCAGCAGCTTCGGCCATAATTTTGTTTTCTAGTATAACTAATGCATGAGAGACAGGCACAT';
dnaseq(24).Header='Anolis_paternus'   ;   dnaseq(24).Sequence='ATGAGCCCATTTATTTTTACAATTATTTTTTCAAGCTTAGCAACAGGCACAATTATTACTATAACCAGCTACCACTGACTTATAGCCTGAGTTGGATTAGAAATAAACACACTAGCAATTATTCCCATTATCTCAACAACACATCACCCACGAGCCACAGAAGCTTCCACAAAATATTTTCTTACACAAGCTGCAGCCTCTGCCTTAATTTTATTTTCTAGTATAACCAATGCATGACATACGGGATCTT';
dnaseq(25).Header='Anolis_sagrei'      ;  dnaseq(25).Sequence='ATGAGCCCAACAATCTTTATAATTATCATACTAAGTCTTGCAACTGGTACAATTATTACTACTACTAGCCACCACTGACTATTAGCCTGAATCGGCCTAGAAATTAATACCCTCTCAATTATTCCAATTATTTCAATAACCCACCACCCACGATCCACAGAAGCCGCTACCAAGTACTTTCTGACACAAGCAGCAGCCTCCGCCCTAATTTTATTTTCCAGTATAACTAATGCATGAGAAACAGGTACAT';
dnaseq(26).Header='Anolis_strahmi'     ;  dnaseq(26).Sequence='ATGAGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACTATAACAAGCTACCACTGACTCATAGCCTGAATTGGACTAGAAATTAACACCTTAGCCATTATTCCAATTATTTCTATACAACACAACCCACGATCAACGGAAGCCGCAACAAAATACTTCCTAACTCAAGCAGCTGCATCCTCCCTAATTTTATTCGCAAGCCTATTTAACGCCTGACAAGTAGGCACAT';
dnaseq(27).Header='Anolis_stratulus'   ;  dnaseq(27).Sequence='ATGAGCCCTATAATTTACACAATCATTTTGTCAAGCCTAGCAACAGGGACAATTATTACCATAACCAGCTACCACTGATTAATAGCTTGAATAGGCCTAGAACTTAATACTCTAGCAATTATTCCATTATCTCATACAACACACAACCCACGATCTACAGAAGCCGCAACAAAATACTTCTTAACACAAGCAGCAGCATCTGCTATAATCTTGTTTTCCAGCATAACCAACGCCTGCTTTACGGGCATAT';
dnaseq(28).Header='Anolis_valencienni' ;  dnaseq(28).Sequence='NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATAACAAGTTATCATTGACTTATGGCCTGAATTGGACTAGAAATAAACACATTAGCAATTATTCCAATTANATTTACTATACACAACCCACGAGCCACAGAAGCTGCCACAAAGTACTTCCTAACACAAGCATCAGCATCAGCCATAATTTTATTTTCAAGCACAATTAACGCCTGACAAACAGGGACCT';
dnaseq(29).Header='Anolis_vanidicus'   ;  dnaseq(29).Sequence='ATGAGCCCAACAATTTATACAACTATTCTAACCAGTCTTGCCACCGGAACAATTATTACAATAACCAGCCACCACTGATTAATAGCCTGANTCGGATTAGAAATNAATACATTAGCAATAATCCCAACTATCTCAACAATGCACCACCCTCGATCAACTGAAGCTGCCACAAAATACTTCTTAACTCAGGCAGCTGCCTCAGCCTTAATTCTCTTTTCTAGTATAACTAACGCCTGACAAACAGGCTCCT';
dnaseq(30).Header='Diplolaemus_darwinii'; dnaseq(30).Sequence='ATGAGCCCAACTACAATAATAATTATTACATCTAGCCTAGCCACGGANACAATCATCACCGCATCAAGCTACCACTGACTACTGGCCTGAGTAGGCCTAGAACTAAATACACTAGCAATTCTTCCAATAATTTCAAAATATCACCACCCACGAGCAACAGAAGCTGCAACAAAATATTTCCTAACACAAGCAGCAGCCTCCGCCATAATCATATTTTCAAGCACACTAAACGCCTGACAAACAGGCACAT';

for ds=1:numel(dnaseq)
    dnaseq(ds).seqvect=nucleo2mat(dnaseq(ds).Sequence) ;
end

% bioinformatics toolbox
SeqsMultiAligned = multialign(dnaseq);
seqalignviewer(SeqsMultiAligned);
1-mean(seqpdist(seqs,'Method','p-distance')); % average pariwise distance
tree = seqlinkage(seqpdist(dnaseq),'single',dnaseq);
phytreeviewer(tree);

% use dtwave to align
addpath('/home/louis/Documents/Matlab/mfiles/DTWaverage') ;
addpath('/home/louis/Documents/Matlab/mfiles/dtw_nucleotide') ;
%[ dist, W ] = DTWaverage(dnaseq(1).seqvect,dnaseq(2).seqvect,[], 0.5, 0, 0) ;
[ dist, W ] = DTWaverage(dnaseq(1).seqvect(:,1:100),dnaseq(2).seqvect(:,20:40),[], 0.5, 0, 1) ;
% get pairwise distances from DTWave
pwdist=zeros(size(dnaseq,2));
for n=1:size(dnaseq,2)
    for m=1:size(dnaseq,2)
        pwdist(n,m) = DTWaverage(dnaseq(n).seqvect,dnaseq(m).seqvect,[], 0.5, 0, 0) ;
        pwdist(m,n) = pwdist(n,m) ;
    end
end
tree = seqlinkage(pwdist,'single',dnaseq);
phytreeviewer(tree);

% create error free 150 long reads (Anolis sequences are 250 long)
reads=[];
for n=1:length(dnaseq)
    for m=1:100
        startpos=randi([1 length(dnaseq(n).Sequence)-150],1,1);
        endpos=startpos+150;
        reads(m+(n-1)*100).seqvect=dnaseq(n).seqvect(:,startpos:endpos);
    end
end



%% generate one 5k long sequence to recover from N reads
addpath('/home/louis/Documents/Matlab/mfiles');
addpath('/home/louis/Documents/Matlab/mfiles/DTWaverage') ;
addpath('/home/louis/Documents/Matlab/mfiles/dtw_nucleotide') ;

haplo(1).Header = 'haplotype1' ;
nucleotides = char(['A' 'C' 'T' 'G']) ;
haplo(1).Sequence = nucleotides(ceil(4*rand(1,5000))) ;
haplo(1).seqvect = nucleo2mat(haplo(1).Sequence) ;

% create error free 150 reads, 1k reads,one 5k sequence = 30x coverage
readnum=1000;
reads=[];
pos=zeros(1,readnum);
for m=1:readnum
    startpos=randi([1 length(haplo(1).Sequence)-150],1,1);
    pos(m)=startpos;
    endpos=startpos+150;
    reads(m).Header=[num2str(m) '_' num2str(startpos)] ;
    reads(m).seqvect=haplo(1).seqvect(:,startpos:endpos);
end

% classification parameters
epoch = 1 ;
LR = [0.9 0.01] ; % learning rate
NS = [1 1] ; % neighboring strength
NC = [1 1] ; % number of children
thexpan = +Inf; % round(0.05*numel(dnaseq)) ;
gama = 0.95 ;
ce = 0 ; % use compression/expansion
fe = length(haplo(1).Sequence) ; % use ends-free alignment
tic; [tree, weight] = ETDTWrec(reads,epoch,LR,NS,NC,thexpan,gama,'seqvect',[],0,ce,fe) ; toc;
%[distSN, distNN] = distSNN(reads,weight,'','seqvect') ;
%depth=5;
%[dico, diconodes] = etreedico(tree,weight,[],[depth 0 0 0 0],[],distSN,distNN,[],1) ;
[score, alignment] = swalign(mat2nucleo(weight{2}),mat2nucleo(haplo(1).seqvect));

% initialise the tree root with reference sequence (needs a reference to guide alignment of reads)
ref = mutatematseq(haplo(1).seqvect,0.3) ; % 70% similarity
tic; [tree, weight] = ETDTWrec(reads,epoch,LR,NS,NC,thexpan,gama,'seqvect',[],0,ce,fe,ref) ; toc;
% visualize with Bioinformatics Toolbox 
sek(1).Sequence=mat2nucleo(ref); sek(1).Header='ref';
sek(2).Sequence=mat2nucleo(haplo(1).seqvect); sek(2).Header='haplo';
sek(3).Sequence=mat2nucleo(weight{2}); sek(3).Header='weight';
seqalignviewer(multialign(sek)) ; % multiple seq alignment
fastawrite('~/Documents/Projects/ShortReads/DTWave_nucleotide/test_one.fa',sek); % to open in clustalx
% WATCH OUT SEEMS TO HAVE BUG WITH INDEX OF SEQUENCES: discrepancy between swalign and seqpdist scores
[score, alignment] = swalign(sek(1),sek(2)) ; showalignment(alignment); % pairwise alignment
seqpdist(sek,'method','p-distance') % pairwise distances
tree = seqlinkage(seqpdist(sek),'single',sek); phytreeviewer(tree);


%% test several haplotypes

% generate haplotypes
nhaplo = 5 ; % 5 haplotypes
len = 5000 ; % genome length 5000bp
haplo=[];
n=1;
haplo(n).Header = ['haplotype' num2str(n)] ;
haplo(n).seqvect = randmatseq(len) ;
haplo(n).Sequence = mat2nucleo(haplo(n).seqvect) ;
for n=2:nhaplo
    haplo(n).Header = ['haplotype' num2str(n)] ;
    haplo(n).seqvect = mutatematseq(haplo(1).seqvect, 0.1) ; % 10% similarity
    haplo(n).Sequence = mat2nucleo(haplo(n).seqvect) ;
end

% generate error free reads
readlen = 150 ; % 150bp long
readnum = (nhaplo * len * 30) / readlen ; % 30x coverage
reads=[];
for m=1:readnum
    haplo_id=randi([1 nhaplo],1,1) ;
    startpos=randi([1 length(haplo(haplo_id).Sequence)-readlen],1,1);
    endpos=startpos+readlen-1;
    reads(m).Header=['read' num2str(m) '_haplo' num2str(haplo_id) '_pos' num2str(startpos)] ;
    reads(m).seqvect=haplo(haplo_id).seqvect(:,startpos:endpos);
end

% dtwave_cluster
reference = mutatematseq(haplo(1).seqvect,0.3) ; % use 30% similarity from first haplotype as reference sequence
epoch = 5 ;
LR = [0.9 0.01] ; % learning rate
NS = [3 3] ; % neighboring strength
NC = [2 2] ; % number of children
thexpan = round(numel(reads)/3) ;
gama = 0.95 ;
ce = 0 ; % use compression/expansion
fe = length(haplo(1).Sequence) ; % use free-ends alignment
tic; [treeA, weightA] = ETDTWrec(reads,1,LR,NS,NC,thexpan,gama,'seqvect',[],0,ce,fe,reference) ; toc;
tic; [treeB, weightB] = ETDTWrec(reads,5,LR,NS,NC,thexpan,gama,'seqvect',[],0,ce,fe,reference) ; toc;
tic; [treeC, weightC] = ETDTWrec(reads,5,[0.01 0.9],NS,NC,thexpan,gama,'seqvect',[],0,ce,fe,reference) ; toc;
tic; [treeD, weightD] = ETDTWrec(reads,5,LR,NS,[4 4],thexpan,gama,'seqvect',[],0,ce,fe,reference) ; toc;
tic; [treeE, weightE] = ETDTWrec(reads,20,LR,NS,[4 4],thexpan,gama,'seqvect',[],0,ce,fe,reference) ; toc;
save('~/Documents/Projects/ShortReads/DTWave_nucleotide/tests_five.dat');

% plot results
% add the weight matrices in the haplotype sequence structure
for w=1:numel(weight)
    haplo(nhaplo+w).Header=['weight' num2str(w)] ;
    haplo(nhaplo+w).seqvect=weight{w};
    haplo(nhaplo+w).Sequence=mat2nucleo(weight{w});
end
% visualize with Bioinformatics Toolbox BUG: RESULTS DIFFERENT THAN CLUSTALX
seqalignviewer(multialign(haplo)) ; % multiple seq alignment
fastawrite('~/Documents/Projects/ShortReads/DTWave_nucleotide/test_five.fa',haplo); % to open in clustalx
% WATCH OUT SEEMS TO HAVE BUG WITH INDEX OF SEQUENCES: discrepancy between swalign and seqpdist scores
phytree = seqlinkage(seqpdist(haplo),'single',haplo); phytreeviewer(phytree);


%% test 2 haplotypes
addpath('/home/louis/Documents/Matlab/mfiles');
addpath('/home/louis/Documents/Matlab/mfiles/DTWaverage') ;
addpath('/home/louis/Documents/Matlab/mfiles/dtw_nucleotide') ;

% generate haplotypes
nhaplo = 2 ; % 2 haplotypes
len = 5000 ; % genome length 5000bp
haplo=[];
n=1;
haplo(n).Header = ['haplotype' num2str(n)] ;
haplo(n).seqvect = randmatseq(len) ;
haplo(n).Sequence = mat2nucleo(haplo(n).seqvect) ;
for n=2:nhaplo
    haplo(n).Header = ['haplotype' num2str(n)] ;
    haplo(n).seqvect = mutatematseq(haplo(1).seqvect, 0.1) ; % 10% similarity
    haplo(n).Sequence = mat2nucleo(haplo(n).seqvect) ;
end

% generate error free reads
readlen = 150 ; % 150bp long
readnum = (nhaplo * len * 30) / readlen ; % 30x coverage
reads=[];
for m=1:readnum
    haplo_id=randi([1 nhaplo],1,1) ;
    startpos=randi([1 length(haplo(haplo_id).Sequence)-readlen],1,1);
    endpos=startpos+readlen-1;
    reads(m).Header=['read' num2str(m) '_haplo' num2str(haplo_id) '_pos' num2str(startpos)] ;
    reads(m).seqvect=haplo(haplo_id).seqvect(:,startpos:endpos);
end

% dtwave_cluster
reference = mutatematseq(haplo(1).seqvect, 0.3) ; % use a third haplotype as a reference
epoch = 5 ;
LR = [0.9 0.01] ; % learning rate
NS = [3 3] ; % neighboring strength
NC = [2 2] ; % number of children
thexpan = round(numel(reads)) ;
gama = 0.95 ;
ce = 0 ; % use compression/expansion
fe = length(haplo(1).Sequence) ; % use free-ends alignment
tic; [treeA, weightA] = ETDTWrec(reads,epoch,LR,NS,NC,thexpan,gama,'seqvect',[],0,ce,fe,reference) ; toc
tic; [treeB, weightB] = ETDTWrec(reads,epoch,LR,NS,NC,2*numel(reads),gama,'seqvect',[],0,ce,fe,reference) ; toc
tic; [treeF, weightF] = ETDTWrec(reads,20,LR,NS,[4 4],round(numel(reads)),gama,'seqvect',[],0,ce,fe,reference) ; toc;
tic; [treeG, weightG] = ETDTWrec(reads,20,LR,NS,[4 4],30*round(numel(reads)),gama,'seqvect',[],0,ce,fe,reference) ; toc;
save('~/Documents/Projects/ShortReads/DTWave_nucleotide/tests_two.dat');

% test ETDTWrec, problem with NS, temporary fix, set NS to single small
% value (0.4), LR to high value not decreasing and use 1-H in the code
% x=0.4; [ exp( (-0^2)/(2*(x)^2) ) exp( (-1^2)/(2*(x)^2) ) exp( (-2^2)/(2*(x)^2) ) exp( (-3^2)/(2*(x)^2) ) ]
% ans = 1.0000    0.0439    0.0000    0.0000
% (for LR=0.7) : 1 - [[1.0000    0.0439    0.0000    0.0000].*0.7]
% ans = 0.3000    0.9692    1.0000    1.0000
% Ideally we don't want the first value to be zero
% e.g.:
tic; [treeH, weightH] = ETDTWrec(reads,2,[0.7 0.7],0.4,[2 2],+Inf,0.95,'seqvect',[],0,0,1,reference) ; toc;

% align one read to reference and parse the position of the read in the
% sequence it was drawn from to compare where it aligns
n=250;
original_pos = str2num(reads(n).Header(strfind(reads(n).Header,'pos')+3:end))+150 ;
[dist, wa, aligned_pos]=DTWaverage(reference,reads(n).seqvect,[], 1, 0, 1) ;
fprintf(1, ['position:' num2str(original_pos) ' - aligned:' num2str(aligned_pos) '\n']);


%% test loading short reads file
addpath('/home/louis/Documents/Matlab/mfiles'); % for ETDTWrec()
addpath('/home/louis/Documents/Matlab/mfiles/DTWaverage') ;
addpath('/home/louis/Documents/Matlab/mfiles/dtw_nucleotide') ;

% read short sequencing reads file
reads=fastqread('/home/louis/Documents/Projects/ShortReads/test_beast_rnd_haplo/1PopDNA_1_9_p1p2.fq');
for m=1:numel(reads)
    reads(m).seqvect=nucleo2mat(reads(m).Sequence);
end

% select and mutate a true haplotype to use as reference
true_seq = fastaread('/home/louis/Documents/Projects/ShortReads/test_beast_rnd_haplo/1PopDNA_1_9.arp.fa.fa');
reference = mutatematseq(nucleo2mat(true_seq(1).Sequence), 0.3) ;

% generate weight sequences
tic; [treeH, weightH] = ETDTWrec(reads,2,[0.7 0.7],0.4,[2 2],numel(reads)/10,0.95,'seqvect',[],0,0,1,reference) ; toc;

% write the tip weights to fasta file
tips = setdiff(1:length(treeH),treeH);
weight_seq=[];
s=1;
for n=tips
    weight_seq(s).Header = ['weight' num2str(n)] ;
    weight_seq(s).seqvect = weightH{n} ;
    weight_seq(s).Sequence = mat2nucleo(weightH{n}) ;
    s=s+1;
end
phytree = seqlinkage(seqpdist([haplo weight_seq]),'average',[haplo weight_seq]); phytreeviewer(phytree); % UPGMA
seqalignviewer(multialign([haplo weight_seq])) ;
fastawrite('/home/louis/Documents/Projects/ShortReads/test_beast_rnd_haplo/1PopDNA_1_9.arp.fa.fa.dtwaven.fa',weight_seq);


%% test 5 short haplotypes
addpath('/home/louis/Documents/Matlab/mfiles');
addpath('/home/louis/Documents/Matlab/mfiles/DTWaverage') ;
addpath('/home/louis/Documents/Matlab/mfiles/dtw_nucleotide') ;

% generate haplotypes
nhaplo = 5 ; % 5 haplotypes
len = 500 ; % genome length 500bp
reference = randmatseq(len) ;
for n=1:nhaplo
    haplo(n).Header = ['haplotype' num2str(n)] ;
    haplo(n).seqvect = mutatematseq(reference, 0.3) ; % 70% dissimilarity
    haplo(n).Sequence = mat2nucleo(haplo(n).seqvect) ;
end

% generate error free reads
readlen = 15 ; % short read length
readnum = (nhaplo * len * 30) / readlen ; % 30x coverage
for m=1:readnum
    haplo_id=randi([1 nhaplo],1,1) ;
    startpos=randi([1 length(haplo(haplo_id).Sequence)-readlen+1],1,1);
    endpos=startpos+readlen-1;
    reads(m).Header=['read' num2str(m) '_haplo' num2str(haplo_id) '_pos' num2str(startpos)] ;
    reads(m).seqvect=haplo(haplo_id).seqvect(:,startpos:endpos);
end
% calculate coverage
coverage=zeros(1,len);
for n=1:readnum
    position=str2double(reads(n).Header(strfind(reads(n).Header,'pos')+3:end)) ; 
    for m=position:(position+readlen-1)
        coverage(m)=coverage(m)+1;
    end
end
bar(coverage);

% generate weight sequences
tic; [treeH, weightH] = ETDTWrec(reads,10,[0.01 0.01],0.4,[2 2],numel(reads)/2,0.95,'seqvect',[],0,0,1,reference) ; toc;

% Visualization
tips = setdiff(1:length(treeH),treeH);
weight_seq=[];
s=1;
for n=tips
    weight_seq(s).Header = ['weight' num2str(n)] ;
    weight_seq(s).seqvect = weightH{n} ;
    weight_seq(s).Sequence = mat2nucleo(weightH{n}) ;
    s=s+1;
end
[~,alignment]=nwalign(haplo(1).Sequence, weight_seq(1).Sequence) ; showalignment(alignment); % pairwise alignment
seqalignviewer(multialign([haplo weight_seq])) ; % multiple seq alignment
phytree = seqlinkage(seqpdist([haplo weight_seq]),'average',[haplo weight_seq]); phytreeviewer(phytree);

% add the weight seqeunces to the "haplo" structure to write them all in the same fasta file
for w=1:numel(weightH)
    haplo(nhaplo+w).Header=['weight' num2str(w)] ;
    haplo(nhaplo+w).seqvect=weightH{w};
    haplo(nhaplo+w).Sequence=mat2nucleo(weightH{w});
end
fastawrite('~/Documents/Projects/ShortReads/DTWave_nucleotide/test_five.fa',haplo); % to open in clustalx


%% testing extra row
reads(2).Header='read1_1del_pos22';
reads(1).seqvect=nucleo2mat('GAGAAATAT');
reads(2).Header='read2_1sub_pos22';
reads(2).seqvect=nucleo2mat('GAGACAATAT');
reference=nucleo2mat('GGCCTCTGTCAATCCGAATCCGAGAAAATATCCTTGTAGTCGCGTAAGCG');
[ dist, updated_weight, aligned_pos ] = DTWaverage(reference, reads(2).seqvect, 1, 0.99, 0, 1 ) ;
weight{1} = [ reference; ones(1,size(reference,2)) ] ;
[ dist, updated_weight, aligned_pos ] = DTWaverage(weight{1}, reads(2).seqvect, 1, 0.99, 0, 1 ) ;
[ dist, updated_weight, aligned_pos ] = DTWaverage(weight{1}, nucleo2mat('GAGAAATATC'), 1, 0.99, 0, 1 ) ;
[ dist, updated_weight, aligned_pos ] = DTWaverage(weight{1}, nucleo2mat('GAGAAAAATATC'), 1, 0.99, 0, 1 ) ;
[ dist, updated_weight, aligned_pos ] = DTWaverage(nucleo2mat('GCATACTCG',1), nucleo2mat('TAACT'), 1, 0.99, 0, 1 ) ; % 2 possibilities if indels same cost as 2 substitutions?
 

%% test subsample etree to estimate N

% read short sequencing reads file
reads=fastqread('/home/louis/Documents/Projects/ShortReads/test_beast_rnd_haplo/1PopDNA_100_p1p2.fq');
for m=1:numel(reads)
    reads(m).seqvect=nucleo2mat(reads(m).Sequence);
end

% select and mutate a true haplotype to use as reference
true_seq = fastaread('/home/louis/Documents/Projects/ShortReads/test_beast_rnd_haplo/1PopDNA_100_1_1.arp.fasta.fa');
reference = mutatematseq(nucleo2mat(true_seq(1).Sequence), 0.3) ;

% generate weight sequences
tic; [treeH, weightH] = ETDTWrec(reads,2,[0.1 0.1],0.4,[2 2],numel(reads)/4,0.95,'seqvect',[],0,0,1,reference) ; toc;
% final etree has only 5 tips
tic; [treeH, weightH] = ETDTWrec(reads,1,[0.1 0.1],0.4,[2 2],numel(reads)/10,0.95,'seqvect',[],0,0,1,reference) ; toc;
% 1816.179284 seconds. final etree has 9 tips (Mapping Precision 0.00408644): 1PopDNA_100_1_1.arp.fasta.fa.dtwaven1.fa
tic; [treeH, weightH] = ETDTWrec(reads,5,[0.1 0.1],0.4,[2 2],numel(reads)/4,0.95,'seqvect',[],0,0,1,reference) ; toc;
% Elapsed time is 10061.380402 seconds. final etree has 9 tips (Mapping Precision 0.00176047): 1PopDNA_100_1_1.arp.fasta.fa.dtwaven2.fa
tic; [treeH, weightH] = ETDTWrec(reads,1,[0.1 0.1],0.4,[2 2],numel(reads)/20,0.95,'seqvect',[],0,0,1,reference) ; toc;
% Elapsed time is 2883.477002 seconds. final etree has 17 tips (Mapping Precision 0.00325232): 1PopDNA_100_1_1.arp.fasta.fa.dtwaven3.fa

% get the weights to run BEAST
tips = setdiff(1:length(treeH),treeH);
weight_seq=[];
s=1;
for n=tips
    weight_seq(s).Header = ['weight' num2str(n)] ;
    weight_seq(s).seqvect = weightH{n} ;
    weight_seq(s).Sequence = mat2nucleo(weightH{n}) ;
    s=s+1;
end
fastawrite('~/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_100_1_1.arp.fasta.fa.dtwaven1.fa',weight_seq);
% clustalw -INFILE=1PopDNA_100_1_1.arp.fasta.fa.dtwaven.fa -OUTPUT=FASTA -OUTFILE=1PopDNA_100_1_1.arp.fasta.fa.dtwaven.fa.aln.fa
% ./beautigen.sh 10 0 1PopDNA_100_1_1.arp.fasta.fa.dtwaven.fa.aln.fa
% for n in {1..100}; do ./beautigen.sh 10 pathsfastafile.aln.fa; done


%% several instances of the same clustering parameters
%
addpath('/home/louis/Documents/Matlab/mfiles');
addpath('/home/louis/Documents/Matlab/mfiles/DTWaverage') ;
addpath('/home/louis/Documents/Matlab/mfiles/dtw_nucleotide') ;

% generate haplotypes
hlen=500; nhaplo=10;
[ haplo, reference ] = sim_haplotypes(nhaplo,hlen,0.3) ;

% number of variable sites
numvarsites(haplo) ;

% generate error free reads
reads = sim_reads(haplo,15,30) ;

% generate weight sequences
tic; [treeH1, weightH1] = ETDTWrec(reads,10,[0.01 0.01],0.4,[2 2],numel(reads),0.95,'seqvect',[],0,0,1,reference) ; toc;
tic; [treeH2, weightH2] = ETDTWrec(reads,10,[0.01 0.01],0.4,[2 2],numel(reads),0.95,'seqvect',[],0,0,1,reference) ; toc;
tic; [treeH3, weightH3] = ETDTWrec(reads,10,[0.01 0.01],0.4,[2 2],numel(reads),0.95,'seqvect',[],0,0,1,reference) ; toc;

L=nhaplo;
for w=1:numel(weightH1)
    haplo(L+w).Header=['weight1_' num2str(w)]; haplo(L+w).seqvect=weightH1{w}; haplo(L+w).Sequence=mat2nucleo(weightH1{w}); end
L=L+numel(weightH1);
for w=1:numel(weightH2)
    haplo(L+w).Header=['weight2_' num2str(w)]; haplo(L+w).seqvect=weightH2{w}; haplo(L+w).Sequence=mat2nucleo(weightH2{w}); end
L=L+numel(weightH2);
for w=1:numel(weightH3)
    haplo(L+w).Header=['weight3_' num2str(w)]; haplo(L+w).seqvect=weightH3{w}; haplo(L+w).Sequence=mat2nucleo(weightH3{w}); end
fastawrite('~/Documents/Projects/ShortReads/DTWave_nucleotide/test_rep.fa',haplo); % to open in clustalx

% plot mapping precision
tic; [tree, weight, BMU, features] = ETDTWrec(reads,30,[0.01 0.01],0.4,[2 2],numel(reads),0.95,'seqvect',[],0,0,1,reference) ; toc;
figure; plot(diff(features.MappingPrecision),'-o') ;
set(gca,'Xtick',1:length(features.NumLeaves),'XTickLabel',features.NumLeaves) ;

%% initialisation of the root with all reads mapped to reference

% load simulated data
% fsc25211 -n 1 -S -i ./1PopDNA_10_5k_short.par
% arp2fasta.awk 1PopDNA_10_5k_short_1_1.arp >1PopDNA_10_5k_short_1_1.fa
% replace_monomorphic ./1PopDNA_10_5k_short_1_1.fa ./1PopDNA_10_5k_short.fasta
% art_illumina -ef -i 1PopDNA_10_5k_short.fasta -l 50 -ss HS25 -f 30 -o single
% cat single_errFree.sam | grep -v ^@ | awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > single_ef.fastq
true_seq = fastaread('/home/louis/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_10_5k_short/1PopDNA_10_5k_short.fasta');
%seqalignviewer(multialign(true_seq)) ;
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
%numvarsites(true_seq) ;
reads=fastqread('/home/louis/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_10_5k_short/single_ef.fastq');
for m=1:numel(reads), reads(m).seqvect=nucleo2mat(reads(m).Sequence) ; end

% generate reference sequence for the root, starting by mutating one of the
reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
% OPTION(1) mutate a true haplotype sequences
reference.seqvect = mutatematseq(true_seq(1).seqvect,0.3) ;
reference.Sequence = mat2nucleo(reference.seqvect) ;
seqpdist([true_seq(1); reference]) ;
numvarsites([true_seq(1); reference]) ; 
% OPTION(2) start with a uniform unknown reference
reference.seqvect = 0.25 + zeros(4,length(true_seq(1).Sequence)) ;

% root weight initialisation: align all reads once to the reference to create initial weight
% choose weight so that, if for a given position all the N reads r are different
% from reference R, then the position in reference is converted to the read
% position value after going through all the reads exactly.
% It cannot reach one so we use 0.9999
% w = (1-0.9999)^(1/N)
w = 0.0001^(1/numel(reads)) ;
rindex = randperm(numel(reads)) ;
reference.seqvect = [ reference.seqvect; ones(1,size(reference.seqvect,2)) ] ; % add persistence vector
for r=rindex % weight must be enough to pull a position toward a different nucleotide 
    [ ~, reference.seqvect ] = DTWaverage( reference.seqvect, reads(r).seqvect, 1, w, 0, 1 ) ;
end
reference.Sequence = mat2nucleo(reference.seqvect) ;
seqalignviewer(multialign([true_seq ; reference])) ; % multiple seq alignment

% generate weight sequences
tic; [tree, weight] = ETDTWrec(reads,40,[1-w 1-w],0.4,[2 2],numel(reads),0.95,'seqvect',[],0,0,1,reference.seqvect(1:4,:)) ; toc;

% Visualization
tips = setdiff(1:length(tree),tree);
weight_seq=[];
s=1;
for n=tips
    weight_seq(s).Header = ['weight' num2str(n)] ;
    weight_seq(s).Sequence = mat2nucleo(weight{n}) ;
    weight_seq(s).seqvect = weight{n} ;
    s=s+1;
end
seqalignviewer(multialign([true_seq ; reference ; weight_seq'])) ; % multiple seq alignment
phytree = seqlinkage(seqpdist([true_seq ; reference ; weight_seq']),'average',[true_seq ; reference ; weight_seq']); phytreeviewer(phytree);
fastawrite('~/Documents/Projects/ShortReads/DTWave_nucleotide/test.fa',[true_seq ; reference ; weight_seq']);

% Results of test with 1PopDNA_10_5k_short dataset, splittig: nbmu>=(father's/2) && entropy>cutoff_95
% 600 sequence matrices loaded
% Parameters -- 10, 0.0152333, 0.0152333, 0.4, 2, 2, 600, 1
% epoch 01/10 21-Jun-2016 12:43:16 -- 0.4, 0.0246116, 0.173763, 5, 3
% epoch 10/10 21-Jun-2016 12:49:01 -- 0.4, 0.00248842, 0.0931555, 637, 319
% Results of test with 1PopDNA_10_5k_short dataset, splittig: nbmu>=father's && entropy>cutoff_95
% epoch 01/10 21-Jun-2016 13:43:33 -- 0.4, 0.0254337, 0.178107, 3, 2
% epoch 10/10 21-Jun-2016 13:44:03 -- 0.4, 0.00654579, 0.0738225, 15, 8
% epoch 01/20 21-Jun-2016 13:46:05 -- 0.4, 0.0254548, 0.177211, 3, 2
% epoch 20/20 21-Jun-2016 13:47:36 -- 0.4, 0.00272706, 0.0590342, 17, 9
% epoch 01/30 21-Jun-2016 13:47:51 -- 0.4, 0.024949, 0.177813, 3, 2
% epoch 30/30 21-Jun-2016 13:50:29 -- 0.4, 0.0020162, 0.0601857, 19, 10
% epoch 01/40 21-Jun-2016 13:51:13 -- 0.4, 0.0246996, 0.175637, 3, 2
% epoch 26/40 21-Jun-2016 13:53:18 -- *** All weights are at minimum entropy95 ***
%  epoch 40/40 21-Jun-2016 15:25:58 -- 0.4, 0.00101375, 0.0585188, 21, 11
% epoch 01/50 21-Jun-2016 14:56:09 -- 0.4, 0.0246132, 0.174728, 3, 2
% epoch 50/50 21-Jun-2016 15:01:09 -- 0.4, 0.00130234, 0.0532278, 21, 11
% epoch 01/75 21-Jun-2016 15:01:16 -- 0.4, 0.0252968, 0.176951, 3, 2
% epoch 75/75 21-Jun-2016 15:09:02 -- 0.4, 0.00406384, 0.0467963, 19, 10
% epoch 001/100 21-Jun-2016 15:09:09 -- 0.4, 0.02541, 0.17894, 3, 2
% epoch 023/100 21-Jun-2016 15:10:53 -- 0.4, 0.00306887, 0.0543982, 13, 7
% epoch 024/100 21-Jun-2016 15:10:58 -- *** All weights are at minimum entropy95 ***
%  epoch 001/100 21-Jun-2016 15:09:09 -- 0.4, 0.02541, 0.17894, 3, 2
%  epoch 100/100 21-Jun-2016 14:14:34 -- 0.4, 0.0029464, 0.0559138, 21, 11
% epoch 001/150 21-Jun-2016 15:27:37 -- 0.4, 0.0245919, 0.174782, 3, 2
% epoch 150/150 21-Jun-2016 15:43:58 -- 0.4, 0.00127617, 0.0477877, 19, 10

% use BEAST to estimate population size
% ../../beautigen.sh 8 0 ./dtwaven10.fa dtwaven10
% ../../beautigen.sh 9 0 ./dtwaven20.fa dtwaven20
% ../../beautigen.sh 8 0 ./dtwaven30.fa dtwaven30
% ../../beautigen.sh 8 0 ./dtwaven40.fa dtwaven40
% ../../beautigen.sh 9 0 ./dtwaven50.fa dtwaven50
% ../../beautigen.sh 10 0 ./dtwaven75.fa dtwaven75
% ../../beautigen.sh 11 0 ./dtwaven100.fa dtwaven100
% ../../beautigen.sh 12 0 ./dtwaven150.fa dtwaven150


%% test for entropy as splitting criterion
addpath('/home/louis/Documents/Matlab/mfiles');
addpath('/home/louis/Documents/Matlab/mfiles/DTWaverage') ;
addpath('/home/louis/Documents/Matlab/mfiles/dtw_nucleotide') ;

true_seq = fastaread('/home/louis/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_10_10k/1PopDNA_10_10k_1_1.arp.fasta.fa');
%seqalignviewer(multialign(true_seq)) ;
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
%numvarsites(true_seq) 
% art_illumina -ef -i 1PopDNA_10_10k_1_1.arp.fasta.fa -l 150 -ss HS25 -f 30 -o single
% cat single_errFree.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > single_ef.fastq;
reads=fastqread('/home/louis/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_10_10k/single_ef.fastq');
for m=1:numel(reads), reads(m).seqvect=nucleo2mat(reads(m).Sequence) ; end

% generate reference sequence for the root, starting by mutating one of the
reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
w = 0.0001^(1/numel(reads)) ;
rindex = randperm(numel(reads)) ;
reference.seqvect = mutatematseq(true_seq(1).seqvect,0.3) ;
reference.seqvect = [ reference.seqvect; ones(1,size(reference.seqvect,2)) ] ; % add persistence vector
for r=rindex % weight must be enough to pull a position toward a different nucleotide 
    [ ~, reference.seqvect ] = DTWaverage( reference.seqvect, reads(r).seqvect, 1, w, 0, 1 ) ;
end
reference.Sequence = mat2nucleo(reference.seqvect) ;

% generate weight sequences
tic; [tree1, weight1] = ETDTWrec(reads,1,[1-w 1-w],0.4,[2 2],numel(reads),1,'seqvect',[],0,0,1,reference.seqvect(1:4,:)) ; toc;
tic; [tree5, weight5] = ETDTWrec(reads,5,[1-w 1-w],0.4,[2 2],numel(reads),1,'seqvect',[],0,0,1,reference.seqvect(1:4,:)) ; toc;
tic; [tree10, weight10] = ETDTWrec(reads,10,[1-w 1-w],0.4,[2 2],numel(reads),1,'seqvect',[],0,0,1,reference.seqvect(1:4,:)) ; toc;
tic; [tree30, weight30] = ETDTWrec(reads,30,[1-w 1-w],0.4,[2 2],numel(reads),1,'seqvect',[],0,0,1,reference.seqvect(1:4,:)) ; toc;
tic; [tree50, weight50] = ETDTWrec(reads,50,[1-w 1-w],0.4,[2 2],numel(reads),1,'seqvect',[],0,0,1,reference.seqvect(1:4,:)) ; toc;
weight2fasta(treeH1,weightH1,'/home/louis/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_10_10k/dtwaven1.fa') ;
weight2fasta(treeH5,weightH5,'/home/louis/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_10_10k/dtwaven5.fa') ;
weight2fasta(treeH10,weightH10,'/home/louis/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_10_10k/dtwaven10.fa') ;
weight2fasta(treeH30,weightH30,'/home/louis/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_10_10k/dtwaven30.fa') ;
weight2fasta(treeH50,weightH50,'/home/louis/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_10_10k/dtwaven50.fa') ;


%% tests with 20 short haplotypes, initialise root, split on shannon entropy, prune tree at each epoch
% cd ~/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_20_5k_short
% fsc25211 -n 1 -S -i ./1PopDNA_20_5k_short.par
% arp2fasta.awk 1PopDNA_20_5k_short_1_1.arp >1PopDNA_20_5k_short_1_1.fa
% replace_monomorphic ./1PopDNA_20_5k_short_1_1.fa ./1PopDNA_20_5k_short.fasta
% art_illumina -ef -i 1PopDNA_20_5k_short.fasta -l 25 -ss HS25 -f 30 -o single
% cat single_errFree.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > single_ef.fastq
cellfun(@(x) addpath(x),{'~/Documents/Matlab/mfiles' '~/Documents/Matlab/mfiles/DTWaverage' '~/Documents/Matlab/mfiles/dtw_nucleotide'});
true_seq = fastaread('/home/louis/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_20_5k_short/1PopDNA_20_5k_short.fasta');
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
reads=fastqread('/home/louis/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_20_5k_short/single_ef.fastq');
for m=1:numel(reads), reads(m).seqvect=nucleo2mat(reads(m).Sequence) ; end
% generate reference sequence for the root, starting by mutating one of the
reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
reference.seqvect = mutatematseq(true_seq(1).seqvect,0.3) ; % OPTION(1) mutate a true haplotype sequences
reference.Sequence = mat2nucleo(reference.seqvect) ;
% root weight initialisation: align all reads once to the reference to create initial weight
w = 0.0001^(1/numel(reads)) ;
rindex = randperm(numel(reads)) ;
reference.seqvect = [ reference.seqvect; ones(1,size(reference.seqvect,2)) ] ; % add persistence vector
for r=rindex % weight must be enough to pull a position toward a different nucleotide 
    [ ~, reference.seqvect ] = DTWaverage( reference.seqvect, reads(r).seqvect, 1, w, 0, 1 ) ;
end
reference.Sequence = mat2nucleo(reference.seqvect) ;
% generate weight sequences
for n=[10 20 30 40 50 75 100]
    tic; ETDTWrec(reads,n,[1-w 1-w],0.4,[2 2],numel(reads),1,'seqvect',[],0,0,1,reference.seqvect(1:4,:)) ; toc;
end
% results are not great as the number of tips increase with the number of epochs...
% maybe the sequences are not long enough, retry with longer sequences
% maybe this is due to the high number of reads?
% maybe this is due to the read length?


%% tests with 20 short haplotypes, initialise root, split on shannon entropy, prune tree at each epoch
% cd ~/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_20_5k_shortB
% fsc25211 -n 1 -S -i ./1PopDNA_20_5k_short.par
% arp2fasta.awk 1PopDNA_20_5k_short/1PopDNA_20_5k_short_1_1.arp >1PopDNA_20_5k_short_1_1.fa
% replace_monomorphic ./1PopDNA_20_5k_short_1_1.fa ./1PopDNA_20_5k_short.fasta
% art_illumina -ef -i 1PopDNA_20_5k_short.fasta -l 75 -ss HS25 -f 30 -o single
% cat single_errFree.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > single_ef.fastq
cellfun(@(x) addpath(x),{'~/Documents/Matlab/mfiles' '~/Documents/Matlab/mfiles/DTWaverage' '~/Documents/Matlab/mfiles/dtw_nucleotide'});
true_seq = fastaread('/home/louis/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_20_5k_shortB/1PopDNA_20_5k_short.fasta');
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
reads=fastqread('/home/louis/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_20_5k_shortB/single_ef.fastq');
for m=1:numel(reads), reads(m).seqvect=nucleo2mat(reads(m).Sequence) ; end
% generate reference sequence for the root, starting by mutating one of the
reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
reference.seqvect = mutatematseq(true_seq(1).seqvect,0.3) ; % OPTION(1) mutate a true haplotype sequences
reference.Sequence = mat2nucleo(reference.seqvect) ;
% root weight initialisation: align all reads once to the reference to create initial weight
w = 0.0001^(1/numel(reads)) ;
rindex = randperm(numel(reads)) ;
reference.seqvect = [ reference.seqvect; ones(1,size(reference.seqvect,2)) ] ; % add persistence vector
for r=rindex % weight must be enough to pull a position toward a different nucleotide 
    [ ~, reference.seqvect ] = DTWaverage( reference.seqvect, reads(r).seqvect, 1, w, 0, 1 ) ;
end
reference.Sequence = mat2nucleo(reference.seqvect) ;
% generate weight sequences
for n=[30 40 50 75 100 150]
    tic; [tree, weight] = ETDTWrec(reads,n,[1-w 1-w],0.4,[2 2],numel(reads),1,'seqvect',[],0,0,1,reference.seqvect(1:4,:)) ; toc;
    weight2fasta(tree,weight,['./dtwaven' num2str(n) '.fa']) ;
end


%% test how the haplotype frequency can be recovered
% cd ~/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_10_5k_short/
cellfun(@(x) addpath(x),{'~/Documents/Matlab/mfiles' '~/Documents/Matlab/mfiles/DTWaverage' '~/Documents/Matlab/mfiles/dtw_nucleotide'});
true_seq = fastaread('~/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_10_5k_short/1PopDNA_10_5k_short.fasta');
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
% write a new fasta file with different frequencies for the true haplotypes
for n=1:numel(true_seq)
    for m=1:n
        fastawrite('./1PopDNA_10_5k_short_varfreq.fa',true_seq(n));
    end
end
% art_illumina -ef -i 1PopDNA_10_5k_short_varfreq.fa -l 50 -ss HS25 -f 30 -o single_varfreq
% cat single_varfreq_errFree.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > single_varfreq_ef.fastq
reads=fastqread('~/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_10_5k_short/single_varfreq_ef.fastq');
for m=1:numel(reads), reads(m).seqvect=nucleo2mat(reads(m).Sequence) ; end
% generate reference sequence for the root, starting by mutating one of the
reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
reference.seqvect = mutatematseq(true_seq(1).seqvect,0.3) ;
reference.Sequence = mat2nucleo(reference.seqvect) ;
% root weight initialisation: align all reads once to the reference to create initial weight
w = 0.0001^(1/numel(reads)) ;
rindex = randperm(numel(reads)) ;
reference.seqvect = [ reference.seqvect; ones(1,size(reference.seqvect,2)) ] ; % add persistence vector
for r=rindex % weight must be enough to pull a position toward a different nucleotide 
    [ ~, reference.seqvect ] = DTWaverage( reference.seqvect, reads(r).seqvect, 1, w, 0, 1 ) ;
end
reference.Sequence = mat2nucleo(reference.seqvect) ;
% generate weight sequences
tic; [tree, weight] = ETDTWrec(reads,25,[1-w 1-w],0.4,[2 2],numel(reads),0.95,'seqvect',[],0,0,1,reference.seqvect(1:4,:)) ; toc;
weight2fasta(tree,weight,['./dtwaven_varfreq.fa']) ;
% align reads to tips and print the counts for each weight
BMU = zeros(numel(reads),2) ;
rposition = zeros(numel(reads),1) ;
for m=1:numel(reads)
    [ treepath, aligned_pos ] = findBMU(reads(m).seqvect, tree, weight, 0, 1) ;
    rposition(m) = aligned_pos ;
    BMU(m,:) = treepath ;
end
save('~/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_10_5k_short/dtwaven_varfreq_BMU.mat','BMU','rposition','tree','weight');
% plotting side by side
bh=bar( [ sort([1 2 3 4 5 6 7 8 9 10 0 0 0 ]'./sum([1 2 3 4 5 6 7 8 9 10]),'descend') sort(arrayfun(@(x) sum(BMU(:,1)==x),unique(BMU(:,1)))./6600,'descend') ] );
set(bh(2),'EdgeAlpha',0,'BarWidth',1,'FaceColor',[0.3 0.5 0.8]);
set(bh(1),'EdgeAlpha',0,'BarWidth',1,'FaceColor',[0.8 0.6 0.4]);
set(gca,'ticklength',[0 0]);
% plot coverage and entropy
tips = setdiff(1:length(tree),tree);
for t=tips
    coverage(weight{t},50,rposition(BMU(:,1)==t),true) ;
    fig = gcf;
    fig.PaperUnits = 'centimeters';
    fig.PaperPosition = [0 0 20 10];
    fig.PaperSize = [20 10] ;
    print(['coverage_varfreq_weight_' num2str(t)],'-dpdf');
end


%% use variance of coverage as a cutoff for splitting
% cd ~/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_10_5k_short/
cellfun(@(x) addpath(x),{'~/Documents/Matlab/mfiles' '~/Documents/Matlab/mfiles/DTWaverage' '~/Documents/Matlab/mfiles/dtw_nucleotide'});
true_seq = fastaread('~/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_10_5k_short/1PopDNA_10_5k_short.fasta');
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
reads=fastqread('~/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_10_5k_short/single_varfreq_ef.fastq');
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
[ coverage, varcov ] = wcoverage(reference.seqvect(1:4,:),50,rpositionr,true) ;
reference.varcov = varcov;
[tree, weight, BMU, features, nbmue, weightvarcov] = ETDTWrec(reads,25,[1-w 1-w],0.4,[2 2],numel(reads),0.95,'seqvect',[],0,0,1,reference) ;
% plot classification tree
figure; treeplot(tree) ; [x,y] = treelayout(tree); text(x,y,num2str(([1:length(tree)]')));
% find BMU for each read and record position
rposition = zeros(numel(reads),1) ;
for m=1:numel(reads)
    [ treepath, aligned_pos ] = findBMU(reads(m).seqvect, tree, weight, 0, 1) ;
    rposition(m) = aligned_pos ;
    BMU(m,:) = treepath ;
end
% plot coverage and entropy
tips = setdiff(1:length(tree),tree);
for t=tips
    [~,wcov] = wcoverage(weight{t},50,rposition(BMU(:,1)==t),true) ;
    fprintf(1,'Node %d: coverage variance: %f\n',t,wcov) ;
    %fig = gcf; % for printing
    fig=figure; % for display
    fig.PaperUnits = 'centimeters';
    fig.PaperPosition = [0 0 20 10];
    fig.PaperSize = [20 10] ;
    %print(['coverage_varfreq_coveragesplit_weight_' num2str(t)],'-dpdf'); % for printing
end
weight2fasta(tree,weight,['./dtwaven_6to2.fa']) ;


%% Test on the alignment position of the short reads
dirm='~/Documents/Matlab/mfiles' ;
cellfun(@(x) addpath(x),{dirm fullfile(dirm,'/DTWaverage') fullfile(dirm,'/dtw_nucleotide')});
truefasta='/home/louis/Documents/Projects/ShortReads/test_popsize_nci/test/2151_20000_5.haplo.fasta';
filefq='/home/louis/Documents/Projects/ShortReads/test_popsize_nci/test/2151_20000_5.combined.fq';
true_seq = fastaread(truefasta);
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
reads=fastqread(filefq);
for m=1:numel(reads), reads(m).seqvect=nucleo2mat(reads(m).Sequence) ; end
reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
reference.seqvect = mutatematseq(true_seq(1).seqvect,0.3) ;
reference.Sequence = mat2nucleo(reference.seqvect) ;
ref_backup=reference;
weight_seq = [] ;
s=1;
g=1;
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
[tree, weight, ~, ~, ~, ~] = ETDTWrec(reads,10,[1-w 1-w],0.9,[2 2],numel(reads),0.95,'seqvect',[],0,0,1,reference) ;


%% Test with last read realignment for each tip in the tree (~reinforcement)
% only the read nucleotide is considered for each locus vector in the reference to compute distance
% 
dirm='~/Documents/Matlab/mfiles' ;
cellfun(@(x) addpath(x),{dirm fullfile(dirm,'/DTWaverage') fullfile(dirm,'/dtw_nucleotide')});
true_seq = fastaread('/home/louis/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_20_5k_shortB/1PopDNA_20_5k_short.fasta');
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
reads=fastqread('/home/louis/Documents/Projects/ShortReads/DTWave_nucleotide/1PopDNA_20_5k_shortB/single_ef.fastq');
for m=1:numel(reads), reads(m).seqvect=nucleo2mat(reads(m).Sequence) ; end
reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
reference.seqvect = mutatematseq(true_seq(1).seqvect,0.3) ;
reference.Sequence = mat2nucleo(reference.seqvect) ;
ref_backup=reference;
weight_seq = [] ;
s=1;
g=1;
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
[tree, weight, BMU, ~, ~, ~] = ETDTWrec(reads,40,[1-w 1-w],0.9,[2 2],numel(reads),0.95,'seqvect',[],0,0,1,reference) ;
plot_dtwnucleo(tree, weight, true_seq, reference) ;
weightf = fortify(5, reads, weight, w, BMU, 0, 1) ;
plot_dtwnucleo(tree, weightf, true_seq, reference) ;


%% test to classify reads by doing multiple etree
% then cluster them based on distance defined as the number of times they fall together 
% generate (short) haplotype dataset
dirm='~/Documents/Matlab/mfiles' ;
cellfun(@(x) addpath(x),{dirm fullfile(dirm,'/DTWaverage') fullfile(dirm,'/dtw_nucleotide')});
hlen=500; nhaplo=10;
reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
[ true_seq, reference.seqvect ] = sim_haplotypes(nhaplo,hlen,0.03) ;
% number of variable sites
if (numvarsites(true_seq)<(hlen*0.1)), error('Too few variable sites') ; end
% show phylo tree
%phytree = seqlinkage(seqpdist(true_seq),'average',true_seq); phytreeviewer(phytree);
% calculate joint entropy based on true haplotypes
%joint_shannonEnt = shannonEntropy_s({true_seq.seqvect}) ;
% generate error free reads
reads = sim_reads(true_seq,20,20) ;
ref_backup=reference;
read_cocount = zeros(numel(reads)) ; % number of co-clustering
read_ald = zeros(numel(reads)) ; % distance between aligned position
% Generate multiple sets of candidate haplotypes
weight_seq = [] ;
s=1;
parfor (g=1:10, 6)
% for g=1:10 %not parallel
    reference=ref_backup;
    % root weight initialisation: align all reads once to the reference to create initial weight
    w = 0.0001^(1/numel(reads)) ;
    rindex = randperm(numel(reads)) ;
    % Subsample?
    % Boost?
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
    reference.varcov = varcov ;
    [reference.entropy, entval] = shannonEntropy(reference.seqvect) ;
    refent_sd = std(entval) ; % standard deviation of the reference entropy vector
    [tree, weight, BMU, ~, ~, ~] = ETDTWrec(reads,10,[1-w 1-w],0.9,[2 2],numel(reads),0.95,'seqvect',[],0,0,1,reference) ;
    [~,~,readpos] = fortify(0,reference.entropy, reads, weight, w, BMU, 0, 1) ;
    % record alignment distance between reads
    read_ald = read_ald + abs(bsxfun(@minus,readpos',readpos)) ;
    % record reads co-classification (how many times each pair of reads are classified together)
    read_cocount = read_cocount + bsxfun(@eq,BMU(:,1)',BMU(:,1)) ;
    % trim the weights by removing high entropy regions
    weight_trim = cellfun(@(x) entrop_trim(x,reference.entropy+2*refent_sd,1), weight, 'UniformOutput', false ) ;
    % align trimmed weights to the true haplotypes
    plot_dtwnucleo(0, weight_trim, true_seq, reference) ;
end
% results for fortify() tests, on 10 replicates of the same data but different order of reads for reference initialisation
% reference.entropy: 0.038662 0.038665 0.038562 0.038659 0.038764 0.038749 0.038689 0.038654 0.038723 0.038759 
% number of iteration to reach reference entropy: 15 18 16 18 15 16 17 16 17 15
save('/home/louis/vecqua/data1000_r20_30x.mat',read_cocount,read_ald) ;
read_dist = 1 - read_cocount./max(read_cocount(:)) ;
%read_dist2 = exp(-(read_cocount./max(read_cocount(:))).^2) ; for i=1:10000, read_dist2(i,i)=0; end;
%read_dist3 = -log(read_cocount./max(read_cocount(:))) ;
figure; hist(squareform(read_dist),100); % not a lot of signal here...
plot(squareform(read_ald),squareform(read_dist),'o') ;set(gca,'Ydir','reverse');
% cluster reads and find number of haplotypes
% MDS plot of reads
[Y] = cmdscale(read_dist) ; % classical MDS (metric, assume euclidean distances)
%[Y] = mdscale(read_dist,2) ; % non classical MDS, much slower...
c = arrayfun(@(x) str2double(cell2mat(regexp(x.Header,'haplo(\d*)_','tokens','once'))),reads) ; % extract haplotypes ID from read header
gscatter(Y(:,1),Y(:,2),c,[],'o') ;% MDS plot colored by haplotypes
% using dp ()
centroids = densipeak(read_dist) ;
% reconstruct haplotype sequences
%...
% plot phylogenetic tree with true haplotypes
%...


%% Test multiple etree with short haplotypes generated by fastsimcoal
cd '/home/louis/Documents/Projects/ShortReads/test_nucleoveq' ;
addpath('/home/louis/Documents/Matlab/mfiles/nucleoveq');
fname='1PopDNA_20_500_5k_07Nov2016_163646';
truefasta=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.haplo.fasta'];
filefq=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.combined.fq'];
true_seq = fastaread(truefasta);
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
reads=fastqread(filefq);
for m=1:numel(reads), reads(m).seqvect=nucleo2mat(reads(m).Sequence) ; end
plot_dtwnucleo(0, [], true_seq) ;
reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
reference.seqvect = mutatematseq(true_seq(1).seqvect,0.02) ;
% OR RANDOM reference.seqvect = randmatseq(length(true_seq(1).Sequence)) ;
reference.Sequence = mat2nucleo(reference.seqvect) ;
ref_backup=reference;
weight_seq = [] ;
s=1;
g=1;
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
plot_dtwnucleo(tree, weight, true_seq, reference, 1) ;


%% Testing entropy value for different haplotype sets
% load and format data
cd '/home/louis/Documents/Projects/ShortReads/test_nucleoveq' ;
addpath('/home/louis/Documents/Matlab/mfiles/nucleoveq');
fname='1PopDNA_10_500_20k_11Nov2016_141255';
fname='1PopDNA_40_500_20k_11Nov2016_220426';
truefasta=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.haplo.fasta'];
filefq=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.combined.fq'];
true_seq = fastaread(truefasta);
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
reads=fastqread(filefq);
read2true=zeros(numel(reads),1); % record which trueseq XX each read comes from, considering Header pattern is '1_XX-n'
for m=1:numel(reads)
    reads(m).seqvect=nucleo2mat(reads(m).Sequence) ;
    read2true(m)=str2double(reads(m).Header(3:strfind(reads(m).Header,'-')-1));
end
plot_dtwnucleo(0, [], true_seq) ;
% align reads to reference
nrep=100 ;
setsize=2:10;
set_ent=zeros(length(setsize),nrep);
ref_ent=zeros(length(setsize),nrep);
ref_entpos=zeros(length(setsize),nrep);
ref_entunik=zeros(length(setsize),nrep);
set_entunik=zeros(length(setsize),nrep);
ref_ent_clust=zeros(length(setsize),nrep);
n=1;
for s=setsize
    parfor repeat=1:nrep
        rndset = randsample(length(true_seq),s) ;
        readset = reads(ismember(read2true,rndset)) ;
        reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
        reference.seqvect = mutatematseq(true_seq(randsample(rndset,1)).seqvect,0.02) ;
        % OR RANDOM reference.seqvect = randmatseq(length(true_seq(1).Sequence)) ;
        reference.Sequence = mat2nucleo(reference.seqvect) ;
        % root weight initialisation: align all reads once to the reference to create initial weight
        w = 0.0001^(1/numel(readset)) ;
        rindex = randperm(numel(readset)) ;
        reference.seqvect = [ reference.seqvect; ones(1,size(reference.seqvect,2)) ] ; % add persistence vector
        BMUr = zeros(numel(readset),2) ;
        rpositionr = zeros(numel(readset),1) ;
        for r=rindex % weight w must be enough to pull a position toward a different nucleotide
            [ dist, reference.seqvect, aligned_pos ] = DTWaverage( reference.seqvect, readset(r).seqvect, 1, w, 0, 1 ) ;
            rpositionr(r) = aligned_pos ;
            BMUr(r,:) = [1 dist] ;
        end
        reference.Sequence = mat2nucleo(reference.seqvect) ;
        [ ~, varcov ] = wcoverage(reference.seqvect(1:4,:),size(readset(1).seqvect,2),rpositionr,false) ;
        reference.varcov = varcov;
        [reference.entropy, shaent] = shannonEntropy(reference.seqvect) ;
        [set_ent(n,repeat), setent] = shannonEntropy_s({true_seq(rndset).seqvect}) ;
        ref_ent(n,repeat) = reference.entropy ;
        ref_entpos(n,repeat) = sum(shaent>0) ;
        ref_entunik(n,repeat) = numel(unique(shaent)) ;
        set_entunik(n,repeat) = numel(unique(setent)) ;
        ref_ent_clust(n,repeat) = numel(densipeak(squareform(pdist(shaent(shaent>0)')))) ;
    end
    fprintf(1,['%d Set entropy: %.3f [%.3f %.3f] - Set ent. unique: %.3f [%.3f %.3f] - ',...
        'Ref entropy: %.3f [%.3f %.3f] - Ref ent. positive: %.3f [%.3f %.3f] - ',...
        'Ref ent. unique: %.3f [%.3f %.3f] - Ref ent. clusters: %.3f [%.3f %.3f]\n'],s,...
        mean(set_ent(n,:)),prctile(set_ent(n,:),5),prctile(set_ent(n,:),95),...
        mean(set_entunik(n,:)),prctile(set_entunik(n,:),5),prctile(set_entunik(n,:),95),...
        mean(ref_ent(n,:)),prctile(ref_ent(n,:),5),prctile(ref_ent(n,:),95),...
        mean(ref_entpos(n,:)),prctile(ref_entpos(n,:),5),prctile(ref_entpos(n,:),95),...
        mean(ref_entunik(n,:)),prctile(ref_entunik(n,:),5),prctile(ref_entunik(n,:),95),...
        mean(ref_ent_clust(n,:)),prctile(ref_ent_clust(n,:),5),prctile(ref_ent_clust(n,:),95) ) ;
    n=n+1;
end
figure
%edges=0:0.05:1.3863;
edges=0:0.01:0.3;
subplot(5,2,1) ; n=1 ; histogram(set_ent(n,:),edges) ; title([num2str(setsize(n)) ' haplotypes']) ;
subplot(5,2,2) ; histogram(ref_ent(n,:),edges) ; title('reads aligned') ;
subplot(5,2,3) ; n=3 ; histogram(set_ent(n,:),edges) ; title([num2str(setsize(n)) ' haplotypes']) ;
subplot(5,2,4) ; histogram(ref_ent(n,:),edges) ; title('reads aligned') ;
subplot(5,2,5) ; n=5 ; histogram(set_ent(n,:),edges) ; title([num2str(setsize(n)) ' haplotypes']) ;
subplot(5,2,6) ; histogram(ref_ent(n,:),edges) ; title('reads aligned') ;
subplot(5,2,7) ; n=7 ; histogram(set_ent(n,:),edges) ; title([num2str(setsize(n)) ' haplotypes']) ;
subplot(5,2,8) ; histogram(ref_ent(n,:),edges) ; title('reads aligned') ;
subplot(5,2,9) ; n=9 ; histogram(set_ent(n,:),edges) ; title([num2str(setsize(n)) ' haplotypes']) ;
subplot(5,2,10) ; histogram(ref_ent(n,:),edges) ; title('reads aligned') ;
% plot number of variable sites and number of unique entropy values
figure
edges=0:30;
edges1=0:3:300;
subplot(5,2,1) ; n=1 ; histogram(set_entunik(n,:),edges) ; title({[num2str(setsize(n)) ' haplotypes'], 'true sequences'}) ;
subplot(5,2,2) ; histogram(ref_entunik(n,:),edges1) ; title({'unique patterns of entropy', 'reads aligned to reference'}) ;
subplot(5,2,3) ; n=3 ; histogram(set_entunik(n,:),edges) ; title({[num2str(setsize(n)) ' haplotypes'], 'true sequences'}) ;
subplot(5,2,4) ; histogram(ref_entunik(n,:),edges1) ; title({'unique patterns of entropy', 'reads aligned to reference'}) ;
subplot(5,2,5) ; n=5 ; histogram(set_entunik(n,:),edges) ; title({[num2str(setsize(n)) ' haplotypes'], 'true sequences'}) ;
subplot(5,2,6) ; histogram(ref_entunik(n,:),edges1) ; title({'unique patterns of entropy', 'reads aligned to reference'}) ;
subplot(5,2,7) ; n=7 ; histogram(set_entunik(n,:),edges) ; title({[num2str(setsize(n)) ' haplotypes'], 'true sequences'}) ;
subplot(5,2,8) ; histogram(ref_entunik(n,:),edges1) ; title({'unique patterns of entropy', 'reads aligned to reference'}) ;
subplot(5,2,9) ; n=9 ; histogram(set_entunik(n,:),edges) ; title({[num2str(setsize(n)) ' haplotypes'], 'true sequences'}) ;
subplot(5,2,10) ; histogram(ref_entunik(n,:),edges1) ; title({'unique patterns of entropy', 'reads aligned to reference'}) ;
save('/home/louis/Documents/Projects/ShortReads/test_nucleoveq/1PopDNA_10_500_20k_11Nov2016_141255_mapping_experiment.mat',...
    'set_ent','set_entunik','ref_ent','ref_entunik','ref_ent_clust');
% plot Entropy value of reference with increasing number of haplotype and arbitrary threshold value
figure; boxplot(ref_ent',2:10,'Color','b','Symbol','o');ylim([0 .6]);hold on;plot([0 10],[.1 .1]);xlabel('Number of haplotypes'); ylabel('Entropy');
figure; boxplot(ref_ent_clust',2:10,'Color','b','Symbol','o');hold on;plot([0 10],[5 5]);xlabel('Number of haplotypes'); ylabel('Entropy Clusters');
% reference entropy is not a very good predictor of the joint entropy of the haplotypes
figure; 
for m=1:length(setsize)
    subplot(3,3,m) ; plot([0 .6],[0 .6],'color',[0.85 0.33 0.1]); hold on; scatter(set_ent(m,:),ref_ent(m,:),'MarkerEdgeColor',[0 0.45 0.74]);
    xlim([0 .6]); ylim([0 .6]); xlabel('true set','FontSize',10); ylabel('reads on reference','FontSize',10);
    title([num2str(setsize(m)) ' haplotypes']);
end

%% test on calculating the probability of a set of haplotype given the reads
% use a constant numebr of true haplotype and map to several fixed number of weights
% load and format data
cd '/home/louis/Documents/Projects/ShortReads/test_nucleoveq' ;
addpath('/home/louis/Documents/Matlab/mfiles/nucleoveq');
%fname='1PopDNA_5_500_5k_11Nov2016_163829';
fname='1PopDNA_40_500_20k_11Nov2016_220426';
truefasta=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.haplo.fasta'];
filefq=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.combined.fq'];
true_seq = fastaread(truefasta);
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
reads=fastqread(filefq);
read2true=zeros(numel(reads),1); % record which trueseq XX each read comes from, considering Header pattern is '1_XX-n'
for m=1:numel(reads)
    reads(m).seqvect=nucleo2mat(reads(m).Sequence) ;
    read2true(m)=str2double(reads(m).Header(3:strfind(reads(m).Header,'-')-1));
end
% plot_dtwnucleo(0, [], true_seq) ;
% T=struct('tree',[],'weight',cell(1,10),'BMU',[]);
% parfor N=1:10
%     [T(N).tree, T(N).weight, T(N).BMU, ~, ~, ~] = ETDTWrec(reads,10,[1-w 1-w],0.9,[N N],numel(reads),0.95,'seqvect',[],0,0,1,reference,0) ;
%     fprintf(1,'P(%d)=%.3f\n',N,p_SHR(reads, T(N).tree, T(N).weight, 0, 1)) ;
% end
nrep=100 ;
setsize=5 ;
numw=1:10 ; % number of weight matrices to project the reads on
pid=fopen('./outfile.txt','w') ;
fprintf(pid,['Haplotype_entropy,Haplotype_entropy_clusters,Short_read_entropy,Short_read_entropy_cluster,',...
    'P_{1}(S|R),P_{1}(H|R),P_{2}(S|R),P_{2}(H|R),P_{3}(S|R),P_{3}(H|R),P_{4}(S|R),P_{4}(H|R),P_{5}(S|R),P_{5}(H|R),',...
    'P_{6}(S|R),P_{6}(H|R),P_{7}(S|R),P_{7}(H|R),P_{8}(S|R),P_{8}(H|R),P_{9}(S|R),P_{9}(H|R),P_{10}(S|R),P_{10}(H|R)\n']);
for repeat=1:nrep
    s=setsize ; % choose randomly "setsize" true haplotype sequences
    rndset = randsample(length(true_seq),s) ;
    readset = reads(ismember(read2true,rndset)) ;
    reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
    reference.seqvect = mutatematseq(true_seq(randsample(rndset,1)).seqvect,0.02) ;
    % OR RANDOM reference.seqvect = randmatseq(length(true_seq(1).Sequence)) ;
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    % root weight initialisation: align all reads once to the reference to create initial weight
    w = 0.0001^(1/numel(readset)) ;
    rindex = randperm(numel(readset)) ;
    reference.seqvect = [ reference.seqvect; ones(1,size(reference.seqvect,2)) ] ; % add persistence vector
    BMUr = zeros(numel(readset),2) ;
    rpositionr = zeros(numel(readset),1) ;
    for r=rindex % weight w must be enough to pull a position toward a different nucleotide
        [ dist, reference.seqvect, aligned_pos ] = DTWaverage( reference.seqvect, readset(r).seqvect, 1, w, 0, 1 ) ;
        rpositionr(r) = aligned_pos ;
        BMUr(r,:) = [1 dist] ;
    end
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    [ ~, varcov ] = wcoverage(reference.seqvect(1:4,:),size(readset(1).seqvect,2),rpositionr,false) ;
    reference.varcov = varcov;
    [reference.entropy, shaent] = shannonEntropy(reference.seqvect) ;
    [set_entropy, setent] = shannonEntropy_s({true_seq(rndset).seqvect}) ;
    fprintf(pid,'%.4f, %d, %.4f, %d, ',set_entropy,numel(densipeak(squareform(pdist(setent(setent>0)')))),...
                                    reference.entropy,numel(densipeak(squareform(pdist(shaent(shaent>0)'))))) ;
    for n=numw
        [tree, weight, BMU, ~, ~, ~] = ETDTWrec(readset,10,[1-w 1-w],0.9,[n n],numel(readset),0.95,'seqvect',[],0,0,1,reference,0) ;
        [ sum_lproba_SgHR, sum_lproba_HgR  ] = p_SHR(readset, tree, weight, 0, 1) ;
        fprintf(pid,'%.4f, %.4f, ', sum_lproba_SgHR, sum_lproba_HgR) ;
    end
    fprintf(pid,'\n');
end
fclose(pid);


%% tests on recovering the true number of haplotypes
cd '/home/louis/Documents/Projects/ShortReads/test_nucleoveq' ;
addpath('/home/louis/Documents/Matlab/mfiles/nucleoveq');
fname='1PopDNA_40_500_20k_11Nov2016_220426';
truefasta=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.haplo.fasta'];
filefq=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.combined.fq'];
true_seq = fastaread(truefasta);
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
reads=fastqread(filefq);
read2true=zeros(numel(reads),1); % record which trueseq XX each read comes from, considering Header pattern is '1_XX-n'
for m=1:numel(reads)
    reads(m).seqvect=nucleo2mat(reads(m).Sequence) ;
    read2true(m)=str2double(reads(m).Header(3:strfind(reads(m).Header,'-')-1));
end
nrep=100 ;
setsize=10 ;
parfor repeat=1:nrep
    pid=fopen(['./outfile_num_haplo10_14Nov2016_' num2str(repeat) '.txt'],'w') ;
    fprintf(pid,'Haplotype_entropy,Haplotype_entropy_clusters,Short_read_entropy,Short_read_entropy_cluster,P(S|R),P(H|R),Num_Weight\n');
    s=setsize ; % choose randomly "setsize" true haplotype sequences
    rndset = randsample(length(true_seq),s) ;
    readset = reads(ismember(read2true,rndset)) ;
    reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
    reference.seqvect = mutatematseq(true_seq(randsample(rndset,1)).seqvect,0.02) ;
    % OR RANDOM reference.seqvect = randmatseq(length(true_seq(1).Sequence)) ;
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    % root weight initialisation: align all reads once to the reference to create initial weight
    w = 0.0001^(1/numel(readset)) ;
    rindex = randperm(numel(readset)) ;
    reference.seqvect = [ reference.seqvect; ones(1,size(reference.seqvect,2)) ] ; % add persistence vector
    BMUr = zeros(numel(readset),2) ;
    rpositionr = zeros(numel(readset),1) ;
    for r=rindex % weight w must be enough to pull a position toward a different nucleotide
        [ dist, reference.seqvect, aligned_pos ] = DTWaverage( reference.seqvect, readset(r).seqvect, 1, w, 0, 1 ) ;
        rpositionr(r) = aligned_pos ;
        BMUr(r,:) = [1 dist] ;
    end
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    [ ~, varcov ] = wcoverage(reference.seqvect(1:4,:),size(readset(1).seqvect,2),rpositionr,false) ;
    reference.varcov = varcov;
    [reference.entropy, shaent] = shannonEntropy(reference.seqvect) ;
    [set_entropy, setent] = shannonEntropy_s({true_seq(rndset).seqvect}) ;
    fprintf(pid,'%.4f, %d, %.4f, %d, ',set_entropy,numel(densipeak(squareform(pdist(setent(setent>0)')))),...
                                    reference.entropy,numel(densipeak(squareform(pdist(shaent(shaent>0)'))))) ;
    [tree, weight, BMU, ~, ~, ~] = ETDTWrec(readset,10,[1-w 1-w],0.9,[2 2],numel(readset),0.95,'seqvect',[],0,0,1,reference,1) ;
    [ sum_lproba_SgHR, sum_lproba_HgR  ] = p_SHR(readset, tree, weight, 0, 1) ;
    fprintf(pid,'%.4f, %.4f, %d\n', sum_lproba_SgHR, sum_lproba_HgR, numel(unique(BMU(:,1)))) ;
    fclose(pid);
end
% load data
d=zeros(100,7);
n=1;
for k = 1:100
    if exist(['./outfile_num_haplo10_14Nov2016_' num2str(k) '.txt'], 'file') == 2
        dat = importdata(['./outfile_num_haplo_14Nov2016_' num2str(k) '.txt'],',',1) ;
        d(k,:) = dat.data ;
        n=n+1;
    end
end
% metadata:
% outfile_num_haplo_13Nov2016_* : 5 true haplotypes, divide at end of each epoch
% outfile_num_haplo_14Nov2016_* : 5 true haplotypes, divide at end of each epoch if Ent>0.1
% outfile_num_haplo10_14Nov2016_* : 10 true haplotypes, divide at end of each epoch if Ent>0.1
figure;scatter(d(d(:,7)>0,1),d(d(:,7)>0,7));
figure;scatter(d(d(:,7)>0,2),d(d(:,7)>0,7));
figure;scatter(d(d(:,7)>0,3),d(d(:,7)>0,7));
figure;scatter(d(d(:,7)>0,4),d(d(:,7)>0,7));
figure;scatter(d(d(:,7)>0,5),d(d(:,7)>0,7));
figure;scatter(d(d(:,7)>0,6),d(d(:,7)>0,7));


%% test stopping criteria
% only divide at the end of each epoch
cd '/home/louis/Documents/Projects/ShortReads/test_nucleoveq/test_numhaplo' ;
addpath('/home/louis/Documents/Matlab/mfiles/nucleoveq');
fname='1PopDNA_40_500_20k_11Nov2016_220426';
truefasta=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.haplo.fasta'];
filefq=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.combined.fq'];
true_seq = fastaread(truefasta);
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
reads=fastqread(filefq);
read2true=zeros(numel(reads),1); % record which trueseq XX each read comes from, considering Header pattern is '1_XX-n'
for m=1:numel(reads)
    reads(m).seqvect=nucleo2mat(reads(m).Sequence) ;
    read2true(m)=str2double(reads(m).Header(3:strfind(reads(m).Header,'-')-1));
end
nrep=100 ;
parfor repeat=1:nrep
    sizes=[5 10 15] ;
    x=randi(3);
    setsize=sizes(x) ;
    s=setsize ; % choose randomly "setsize" true haplotype sequences
    rndset = randsample(length(true_seq),s) ;
    readset = reads(ismember(read2true,rndset)) ;
    reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
    reference.seqvect = mutatematseq(true_seq(randsample(rndset,1)).seqvect,0.02) ;
    % OR RANDOM reference.seqvect = randmatseq(length(true_seq(1).Sequence)) ;
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    % root weight initialisation: align all reads once to the reference to create initial weight
    w = 0.0001^(1/numel(readset)) ;
    rindex = randperm(numel(readset)) ;
    reference.seqvect = [ reference.seqvect; ones(1,size(reference.seqvect,2)) ] ; % add persistence vector
    BMUr = zeros(numel(readset),2) ;
    rpositionr = zeros(numel(readset),1) ;
    for r=rindex % weight w must be enough to pull a position toward a different nucleotide
        [ dist, reference.seqvect, aligned_pos ] = DTWaverage( reference.seqvect, readset(r).seqvect, 1, w, 0, 1 ) ;
        rpositionr(r) = aligned_pos ;
        BMUr(r,:) = [1 dist] ;
    end
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    [ ~, varcov ] = wcoverage(reference.seqvect(1:4,:),size(readset(1).seqvect,2),rpositionr,false) ;
    reference.varcov = varcov;
    [reference.entropy, shaent] = shannonEntropy(reference.seqvect) ;
    [set_entropy, setent] = shannonEntropy_s({true_seq(rndset).seqvect}) ;
    [tree, weight, BMU, ~, ~, ~] = ETDTWrec(readset,10,[1-w 1-w],0.9,[2 2],numel(readset),0.95,'seqvect',[],0,0,1,reference,1) ;
    AllEntropy = shannonEntropy_s(weight(unique(BMU(:,1))')) ;
    [ sum_lproba_SgHR, sum_lproba_HgR  ] = p_SHR(readset, tree, weight, 0, 1) ;
    pid=fopen(['./outfile_num10_18Nov2016_' num2str(repeat) '.txt'],'w') ;
    fprintf(pid,'%d, %.4f, %.4f, %.4f, %.4f, %d\n',setsize, set_entropy, AllEntropy, sum_lproba_SgHR, sum_lproba_HgR, numel(unique(BMU(:,1)))) ;
    fclose(pid);
end
% load data
d=zeros(100,6);
n=1;
for k = 1:100
    if exist(['./outfile_num10_18Nov2016_' num2str(k) '.txt'], 'file') == 2
        dat = importdata(['./outfile_num10_18Nov2016_' num2str(k) '.txt'],',') ;
        d(k,:) = dat ;
        n=n+1;
    end
end
figure; scatter(d(d(:,1)>0,1),d(d(:,1)>0,2));xlim([0 20]); ylabel('True haplotype set entropy');
figure; scatter(d(d(:,1)>0,1),d(d(:,1)>0,3));xlim([0 20]); ylabel('Reconstructed haplotype set entropy');
figure; scatter(d(d(:,1)>0,5),d(d(:,1)>0,4),d(d(:,1)>0,6).^3,d(d(:,1)>0,1),'filled');xlabel('log(P(H|R))');ylabel('log(P(S|H))');
figure; scatter(d(d(:,1)>0,2),d(d(:,1)>0,3),[],d(d(:,1)>0,1),'filled'); xlabel('True set entropy'); ylabel('Reconstructed set entropy');


%% test 3 haplotypes no dividing reads taken in order
cd '/home/louis/Documents/Projects/ShortReads/test_nucleoveq' ;
addpath('/home/louis/Documents/Matlab/mfiles/nucleoveq');
fname='1PopDNA_40_500_20k_11Nov2016_220426';
truefasta=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.haplo.fasta'];
filefq=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.combined.fq'];
true_seq = fastaread(truefasta);
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
reads=fastqread(filefq);
read2true=zeros(numel(reads),1); % record which trueseq XX each read comes from, considering Header pattern is '1_XX-n'
for m=1:numel(reads)
    reads(m).seqvect=nucleo2mat(reads(m).Sequence) ;
    read2true(m)=str2double(reads(m).Header(3:strfind(reads(m).Header,'-')-1));
end
% T=struct('tree',[],'weight',cell(1,10),'BMU',[]);
% parfor N=1:10
%     [T(N).tree, T(N).weight, T(N).BMU, ~, ~, ~] = ETDTWrec(reads,10,[1-w 1-w],0.9,[N N],numel(reads),0.95,'seqvect',[],0,0,1,reference,0) ;
%     fprintf(1,'P(%d)=%.3f\n',N,p_SHR(reads, T(N).tree, T(N).weight, 0, 1)) ;
% end
nrep=100 ;
setsize=3 ;
numw=3 ; % number of weight matrices to project the reads on
pid=fopen('./pooling3/outfile_pooling_rshuff.txt','w') ;
fprintf(pid,'Pct_sim(weight;best_true),num_SNPs,Pct_sim(true;true)\n');
for repeat=1:nrep
    rndset = randsample(length(true_seq),setsize) ; % choose randomly "setsize" true haplotype sequences
    readset = reads(ismember(read2true,rndset)) ;
    reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
    reference.seqvect = mutatematseq(true_seq(randsample(rndset,1)).seqvect,0.02) ;
    %reference.seqvect = randmatseq(length(true_seq(1).Sequence)) ; % RANDOM
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    % root weight initialisation: align all reads once to the reference to create initial weight
    w = 0.0001^(1/numel(readset)) ;
    rindex = randperm(numel(readset)) ;
    reference.seqvect = [ reference.seqvect; ones(1,size(reference.seqvect,2)) ] ; % add persistence vector
    BMUr = zeros(numel(readset),2) ;
    rpositionr = zeros(numel(readset),1) ;
    for r=rindex % weight w must be enough to pull a position toward a different nucleotide
        [ dist, reference.seqvect, aligned_pos ] = DTWaverage( reference.seqvect, readset(r).seqvect, 1, w, 0, 1 ) ;
        rpositionr(r) = aligned_pos ;
        BMUr(r,:) = [1 dist] ;
    end
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    [ ~, varcov ] = wcoverage(reference.seqvect(1:4,:),size(readset(1).seqvect,2),rpositionr,false) ;
    reference.varcov = varcov;
    [reference.entropy, shaent] = shannonEntropy(reference.seqvect) ;
    [set_entropy, setent] = shannonEntropy_s({true_seq(rndset).seqvect}) ;
    %for n=numw
        [~,I]=sort(rpositionr); % choose the read in order to their position in reference rather than random
        %I=0; % to disable ordering of the reads
        [tree, weight, BMU, ~, ~, ~] = ETDTWrec(readset,1,[1-(w/3) 1-(w/3)],0.01,[numw numw],numel(readset),0.95,'seqvect',[],0,0,1,reference,0,I') ;
        %[ sum_lproba_SgHR, sum_lproba_HgR  ] = p_SHR(readset, tree, weight, 0, 1) ;
        %fprintf(pid,'%.4f, %.4f, ', sum_lproba_SgHR, sum_lproba_HgR) ;
    %end
    %plot_dtwnucleo(tree, weight, true_seq(rndset)) ;
    %pause;
    aveminidist = avedist(tree, weight, true_seq(rndset)) ; % average distance between a weight and its closest true sequence
    Profile = seqprofile(true_seq(rndset),'Alphabet','NT');         % number of SNPs in the true sequences alignment
    fprintf( pid,'%.4f, %d, %.4f\n', aveminidist, length(Profile)-sum(sum(Profile==1)),...
        mean(1-seqpdist(true_seq(rndset),'Method','p-distance')) ) ;
end
fclose(pid);

%% test 3 haplotypes with different frequencies no dividing reads taken in order
cd '/home/louis/Documents/Projects/ShortReads/test_nucleoveq/pooling3' ;
addpath('/home/louis/Documents/Matlab/mfiles/nucleoveq');
fname='1PopDNA_40_500_20k_11Nov2016_220426';
truefasta=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.haplo.fasta'];
true_seq = fastaread(truefasta);
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
nrep=100 ;
setsize=3 ;
numw=3 ; % number of weight matrices to project the reads on
pid=fopen('./outfile_pooling_ralignorder_difffreq.txt','w') ;
fprintf(pid,'Pct_sim(weight;best_true),num_SNPs,Pct_sim(true;true)\n');
for repeat=1:nrep
    rndset = randsample(length(true_seq),setsize) ; % choose randomly "setsize" true haplotype sequences
    % write a new fasta file with different frequencies for the true haplotypes (1:2:5)
    warnState = warning; %Save the current warning state
    warning('off','Bioinfo:fastawrite:AppendToFile');
    fastawrite(['./' fname '_' num2str(repeat) '_varfreq.fa'],true_seq(rndset(1)));
    fastawrite(['./' fname '_' num2str(repeat) '_varfreq.fa'],true_seq(rndset(2)));
    fastawrite(['./' fname '_' num2str(repeat) '_varfreq.fa'],true_seq(rndset(2)));
    fastawrite(['./' fname '_' num2str(repeat) '_varfreq.fa'],true_seq(rndset(3)));
    fastawrite(['./' fname '_' num2str(repeat) '_varfreq.fa'],true_seq(rndset(3)));
    fastawrite(['./' fname '_' num2str(repeat) '_varfreq.fa'],true_seq(rndset(3)));
    fastawrite(['./' fname '_' num2str(repeat) '_varfreq.fa'],true_seq(rndset(3)));
    fastawrite(['./' fname '_' num2str(repeat) '_varfreq.fa'],true_seq(rndset(3)));
    warning(warnState) %Reset warning state to previous settings
    system(['/home/louis/Downloads/art_bin_ChocolateCherryCake/art_illumina -ef -i ./'...
        fname '_' num2str(repeat) '_varfreq.fa -l 50 -ss HS25 -f 30 -o single_varfreq_' num2str(repeat)]) ;
    system(['cat single_varfreq_' num2str(repeat)...
        '_errFree.sam | grep -v ^@ | awk ''{print "@"$1"\n"$10"\n+\n"$11}'' > single_varfreq_' num2str(repeat) '_ef.fastq']) ;
    readset = fastqread(['single_varfreq_' num2str(repeat) '_ef.fastq']);
    read2true=zeros(numel(readset),1); % record which trueseq XX each read comes from, considering Header pattern is '1_XX-n'
    for m=1:numel(readset)
        readset(m).seqvect=nucleo2mat(readset(m).Sequence) ;
        read2true(m)=str2double(readset(m).Header(3:strfind(readset(m).Header,'-')-1));
    end
    reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
    reference.seqvect = mutatematseq(true_seq(randsample(rndset,1)).seqvect,0.02) ;
    %reference.seqvect = randmatseq(length(true_seq(1).Sequence)) ; % RANDOM
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    % root weight initialisation: align all reads once to the reference to create initial weight
    w = 0.0001^(1/numel(readset)) ;
    rindex = randperm(numel(readset)) ;
    reference.seqvect = [ reference.seqvect; ones(1,size(reference.seqvect,2)) ] ; % add persistence vector
    BMUr = zeros(numel(readset),2) ;
    rpositionr = zeros(numel(readset),1) ;
    for r=rindex % weight w must be enough to pull a position toward a different nucleotide
        [ dist, reference.seqvect, aligned_pos ] = DTWaverage( reference.seqvect, readset(r).seqvect, 1, w, 0, 1 ) ;
        rpositionr(r) = aligned_pos ;
        BMUr(r,:) = [1 dist] ;
    end
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    [ ~, varcov ] = wcoverage(reference.seqvect(1:4,:),size(readset(1).seqvect,2),rpositionr,false) ;
    reference.varcov = varcov;
    [reference.entropy, shaent] = shannonEntropy(reference.seqvect) ;
    [set_entropy, setent] = shannonEntropy_s({true_seq(rndset).seqvect}) ;
    %for n=numw
        [~,I]=sort(rpositionr); % choose the read in order to their position in reference rather than random
        %I=0; % to disable ordering of the reads
        [tree, weight, BMU, ~, ~, ~] = ETDTWrec(readset,1,[1-(w/3) 1-(w/3)],0.01,[numw numw],numel(readset),0.95,'seqvect',[],0,0,1,reference,0,I') ;
        %[ sum_lproba_SgHR, sum_lproba_HgR  ] = p_SHR(readset, tree, weight, 0, 1) ;
        %fprintf(pid,'%.4f, %.4f, ', sum_lproba_SgHR, sum_lproba_HgR) ;
    %end
    %plot_dtwnucleo(tree, weight, true_seq(rndset)) ;
    %pause;
    aveminidist = avedist(tree, weight, true_seq(rndset)) ; % average distance between a weight and its closest true sequence
    Profile = seqprofile(true_seq(rndset),'Alphabet','NT'); % number of SNPs in the true sequences alignment
    fprintf( pid,'%.4f, %d, %.4f\n', aveminidist, length(Profile)-sum(sum(Profile==1)),...
        mean(1-seqpdist(true_seq(rndset),'Method','p-distance')) ) ;
end
fclose(pid);



%% test 3 haplotypes no dividing reads taken in order, reads taken with errors
cd '/home/louis/Documents/Projects/ShortReads/test_nucleoveq/pooling3' ;
addpath('/home/louis/Documents/Matlab/mfiles/nucleoveq');
fname='1PopDNA_40_500_20k_11Nov2016_220426';
truefasta=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.haplo.fasta'];
true_seq = fastaread(truefasta);
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
nrep=100 ;
setsize=3 ;
numw=3 ; % number of weight matrices to project the reads on
pid=fopen('./outfile_pooling_ralignorder_errorreads_trimmed.txt','w') ;
fprintf(pid,'Pct_sim(weight;best_true),num_SNPs,Pct_sim(true;true)\n');
for repeat=1:nrep
    rndset = randsample(length(true_seq),setsize) ; % choose randomly "setsize" true haplotype sequences
    warnState = warning; %Save the current warning state
    warning('off','Bioinfo:fastawrite:AppendToFile');
    fastawrite(['./' fname '_' num2str(repeat) '_readerr.fa'],true_seq(rndset(1)));
    fastawrite(['./' fname '_' num2str(repeat) '_readerr.fa'],true_seq(rndset(2)));
    fastawrite(['./' fname '_' num2str(repeat) '_readerr.fa'],true_seq(rndset(3)));
    warning(warnState) %Reset warning state to previous settings
    system(['/home/louis/Downloads/art_bin_ChocolateCherryCake/art_illumina -ef -i ./' fname '_' num2str(repeat) ...
        '_readerr.fa -l 50 -ss HS25 -f 30 -o single_readerr_' num2str(repeat)]) ;
    delete(['./' fname '_' num2str(repeat) '_readerr.fa']);
    % Read quality trimming and filtering
    filtrd = seqtrim(['single_readerr_' num2str(repeat) '.fq'],'Method','MaxNumberLowQualityBases','Threshold',[0 15],'WindowSize',10) ;
    filtrd2 = seqfilter(filtrd,'Method','MinLength','Threshold',1) ;
    readset = fastqread(filtrd2{1});
    %readset = fastqread(['single_readerr_' num2str(repeat) '.fq']);
    for m=1:numel(readset)
        quali = qual2accu(readset(m).Quality,33) ;
        readset(m).seqvect=nucleo2mat(readset(m).Sequence,quali) ;
    end
    reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
    reference.seqvect = mutatematseq(true_seq(randsample(rndset,1)).seqvect,0.02) ;
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    % root weight initialisation: align all reads once to the reference to create initial weight
    w = 0.0001^(1/numel(readset)) ;
    rindex = randperm(numel(readset)) ;
    reference.seqvect = [ reference.seqvect; ones(1,size(reference.seqvect,2)) ] ; % add persistence vector
    BMUr = zeros(numel(readset),2) ;
    rpositionr = zeros(numel(readset),1) ;
    for r=rindex % weight w must be enough to pull a position toward a different nucleotide
        [ dist, reference.seqvect, aligned_pos ] = DTWaverage( reference.seqvect, readset(r).seqvect, 1, w, 0, 1 ) ;
        rpositionr(r) = aligned_pos ;
        BMUr(r,:) = [1 dist] ;
    end
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    %[ ~, varcov ] = wcoverage(reference.seqvect(1:4,:),size(readset(1).seqvect,2),rpositionr,false) ;
    %reference.varcov = varcov;
    reference.varcov = 0;
    [reference.entropy, shaent] = shannonEntropy(reference.seqvect) ;
    [set_entropy, setent] = shannonEntropy_s({true_seq(rndset).seqvect}) ;
    %for n=numw
        [~,I]=sort(rpositionr); % choose the read in order to their position in reference rather than random
        %I=0; % to disable ordering of the reads
        [tree, weight, BMU, ~, ~, ~] = ETDTWrec(readset,1,[1-(w/3) 1-(w/3)],0.01,[numw numw],numel(readset),0.95,'seqvect',[],0,0,1,reference,0,I') ;
        %[ sum_lproba_SgHR, sum_lproba_HgR  ] = p_SHR(readset, tree, weight, 0, 1) ;
        %fprintf(pid,'%.4f, %.4f, ', sum_lproba_SgHR, sum_lproba_HgR) ;
    %end
    %plot_dtwnucleo(tree, weight, true_seq(rndset)) ;
    %pause;
    aveminidist = avedist(tree, weight, true_seq(rndset)) ; % average distance between a weight and its closest true sequence
    Profile = seqprofile(true_seq(rndset),'Alphabet','NT'); % number of SNPs in the true sequences alignment
    fprintf( pid,'%.4f, %d, %.4f\n', aveminidist, length(Profile)-sum(sum(Profile==1)),...
        mean(1-seqpdist(true_seq(rndset),'Method','p-distance')) ) ;
end
fclose(pid);


%% test estimating the original number of haplotypes using reference matrix only
cd '/home/louis/Documents/Projects/ShortReads/test_nucleoveq/pooling3' ;
addpath('/home/louis/Documents/Matlab/mfiles/nucleoveq');
fname='1PopDNA_40_500_20k_11Nov2016_220426';
truefasta=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.haplo.fasta'];
true_seq = fastaread(truefasta);
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
filefq=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.combined.fq'];
reads=fastqread(filefq);
read2true=zeros(numel(reads),1); % record which trueseq XX each read comes from, considering Header pattern is '1_XX-n'
for m=1:numel(reads)
    reads(m).seqvect=nucleo2mat(reads(m).Sequence) ;
    read2true(m)=str2double(reads(m).Header(3:strfind(reads(m).Header,'-')-1));
end
nrep=100 ;
pid=fopen('./outfile_pooling_numhaplo.txt','w') ;
fprintf(pid,'Num_haplo,Estimate\n');
for repeat=1:nrep
    setsize=randi([2 10]) ;
    rndset = randsample(length(true_seq),setsize) ; % choose randomly "setsize" true haplotype sequences
    readset = reads(ismember(read2true,rndset)) ;
    reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
    reference.seqvect = mutatematseq(true_seq(randsample(rndset,1)).seqvect,0.02) ;
    %reference.seqvect = randmatseq(length(true_seq(1).Sequence)) ; % RANDOM
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    % root weight initialisation: align all reads once to the reference to create initial weight
    w = 0.0001^(1/numel(readset)) ;
    rindex = randperm(numel(readset)) ;
    reference.seqvect = [ reference.seqvect; ones(1,size(reference.seqvect,2)) ] ; % add persistence vector
    BMUr = zeros(numel(readset),2) ;
    rpositionr = zeros(numel(readset),1) ;
    for r=rindex % weight w must be enough to pull a position toward a different nucleotide
        [ dist, reference.seqvect, aligned_pos ] = DTWaverage( reference.seqvect, readset(r).seqvect, 1, w, 0, 1 ) ;
        rpositionr(r) = aligned_pos ;
        BMUr(r,:) = [1 dist] ;
    end
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    [ ~, varcov ] = wcoverage(reference.seqvect(1:4,:),size(readset(1).seqvect,2),rpositionr,false) ;
    reference.varcov = varcov;
    densipeak(squareform(pdist(sort(reference.seqvect(1:4,:))')))
    fprintf( pid,'%d, %d\n', setsize, ) ;
end
fclose(pid);
%%% DOES NOT WORK TOO MANY CENTROIDS ARE FOUND


%% test 3/5/10 haplotypes no dividing reads taken in order, longer haplotype sequence and read length
% art_illumina -ef -i 1PopDNA_40_1200_20k_28Nov2016_171414.haplo.fasta -p -l 150 -ss HS25 -f 30 -m 530 -s 0 -o 1PopDNA_40_1200_20k_28Nov2016_171414.paired
cd '/home/louis/Documents/Projects/ShortReads/test_nucleoveq/pooling3' ;
addpath('/home/louis/Documents/Matlab/mfiles/nucleoveq');
fname='1PopDNA_40_1200_20k_28Nov2016_171414';
truefasta=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.haplo.fasta'];
filefq=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.combined.fq'];
true_seq = fastaread(truefasta);
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
reads=fastqread(filefq);
read2true=zeros(numel(reads),1); % record which trueseq XX each read comes from, considering Header pattern is '1_XX-n'
for m=1:numel(reads)
    reads(m).seqvect=nucleo2mat(reads(m).Sequence) ;
    read2true(m)=str2double(reads(m).Header(3:strfind(reads(m).Header,'-')-1));
end
nrep=100 ;
setsize=3 ;
numw=setsize ; % number of weight matrices to project the reads on
%pid=fopen('./outfile_pooling10_ralignorder_long.txt','w') ;
pid=fopen('./outfile_pooling3_ralignorder_long_pctTrueSibling.txt','w') ;
%pid=1;
fprintf(pid,'Pct_sim(weight;best_true),num_SNPs,Pct_sim(true;true),Pct_sim(true_sibling)\n');
for repeat=1:nrep
    rndset = randsample(length(true_seq),setsize) ; % choose randomly "setsize" true haplotype sequences
    readset = reads(ismember(read2true,rndset)) ;
    reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
    reference.seqvect = mutatematseq(true_seq(randsample(rndset,1)).seqvect,0.02) ;
    %reference.seqvect = randmatseq(length(true_seq(1).Sequence)) ; % RANDOM REFERENCE
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    % root weight initialisation: align all reads once to the reference to create initial weight
    w = 0.0001^(1/numel(readset)) ;
    rindex = randperm(numel(readset)) ;
    reference.seqvect = [ reference.seqvect; ones(1,size(reference.seqvect,2)) ] ; % add persistence vector
    BMUr = zeros(numel(readset),2) ;
    rpositionr = zeros(numel(readset),1) ;
    for r=rindex % weight w must be enough to pull a position toward a different nucleotide
        [ dist, reference.seqvect, aligned_pos ] = DTWaverage( reference.seqvect, readset(r).seqvect, 1, w, 0, 1 ) ;
        rpositionr(r) = aligned_pos ;
        BMUr(r,:) = [1 dist] ;
    end
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    [ ~, reference.varcov ] = wcoverage(reference.seqvect(1:4,:),size(readset(1).seqvect,2),rpositionr,false) ;
    [reference.entropy, shaent] = shannonEntropy(reference.seqvect) ;
    %[set_entropy, setent] = shannonEntropy_s({true_seq(rndset).seqvect}) ;
    [~,I]=sort(rpositionr); % choose the read in order to their position in reference rather than random
    %[tree, weight, BMU, ~, ~, ~] = ETDTWrec(readset,1,[1-(w/3) 1-(w/3)],0.01,[numw numw],numel(readset),0.95,'seqvect',[],0,0,1,reference,0,I') ;
    [tree, weight, BMU, ~, ~, ~] = ETDTWrec(readset,1,[1-(w/numw) 1-(w/numw)],0.01,[numw numw],numel(readset),0.95,'seqvect',[],0,0,1,reference,0,I') ;
    [aveminidist, pct_mini_is_true] = avedist(tree, weight, true_seq(rndset)) ; % average distance between a weight and its closest true sequence
    Profile = seqprofile(true_seq(rndset),'Alphabet','NT'); % number of SNPs in the true sequences alignment
    fprintf( pid,'%.4f, %d, %.4f, %.4f\n', aveminidist, length(Profile)-sum(sum(Profile==1)),...
        mean(1-seqpdist(true_seq(rndset),'Method','p-distance')), pct_mini_is_true ) ;
end
fclose(pid);


%% test indels on 1 haplotype
cd '/home/louis/Documents/Projects/ShortReads/test_nucleoveq/pooling3' ;
addpath('/home/louis/Documents/Matlab/mfiles/nucleoveq');
fname='1PopDNA_40_1200_20k_28Nov2016_171414';
truefasta=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.haplo.fasta'];
filefq=['/home/louis/Documents/Projects/ShortReads/test_nucleoveq/' fname '.combined.fq'];
true_seq = fastaread(truefasta);
for m=1:numel(true_seq), true_seq(m).seqvect=nucleo2mat(true_seq(m).Sequence) ; end
reads=fastqread(filefq);
read2true=zeros(numel(reads),1); % record which trueseq XX each read comes from, considering Header pattern is '1_XX-n'
for m=1:numel(reads)
    reads(m).seqvect=nucleo2mat(reads(m).Sequence) ;
    read2true(m)=str2double(reads(m).Header(3:strfind(reads(m).Header,'-')-1));
end
nrep=100 ;
setsize=1 ;
pid=fopen('./outfile_pooling_long_indels.txt','w') ;
fprintf(pid,'Pct_sim(weight;true),length(reference),length(reconstructed),length(true)\n');
for repeat=1:nrep
    rndset = randsample(length(true_seq),setsize) ;
    readset = reads(ismember(read2true,rndset)) ;
    reference = struct('Header','root reference','Sequence','','seqvect',[]) ;
    reference.seqvect = mutatematseq(true_seq(rndset).seqvect,0.1) ;
    % introduce indels in reference
    % In human, 0.04 indels per 100 nucleotides (1.5 million in total in human genome) from:
    % Mills RE, Luttig CT, Larkins CE, Beauchamp A, Tsui C, Pittard WS, Devine SE: An initial map of insertion and deletion (INDEL) variation in the human genome. Genome Research. 2006, 16 (9): 1182-1190. 10.1101/gr.4565806.
    x=[ 1 0 0 0 ]' ; % 1 nucleotide
    numin = randi(5) ; % up to 5 insertion
    inpos = sort( randi(length(reference.seqvect),1,numin) ) ; % choose numin random positions
    for n = 1:numin
        reference.seqvect = [ reference.seqvect(:,1:inpos(n)) [x(randperm(4))] reference.seqvect(:,(inpos(n)+1):end)] ;
    end
    numdel = randi(5) ; % up to 5 deletion
    delpos = sort( randi(length(reference.seqvect)-numdel,1,numdel) ) ;
    for n = 1:numdel
        reference.seqvect = [ reference.seqvect(:,1:(delpos(n)-1)) reference.seqvect(:,(delpos(n)+1):end) ] ;
    end
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    L = length(reference.Sequence) ;
    % align all reads once to the reference to reconstruct haplotype
    w = 0.0001^(1/numel(readset)) ;
    rindex = randperm(numel(readset)) ;
    reference.seqvect = [ reference.seqvect; ones(1,size(reference.seqvect,2)) ] ; % add persistence vector
    for r=rindex % weight w must be enough to pull a position toward a different nucleotide
        [ dist, reference.seqvect, aligned_pos ] = DTWaverage( reference.seqvect, readset(r).seqvect, 1, w, 0, 1 ) ;
        reference.seqvect = reference.seqvect(:, reference.seqvect(5,:)>0.9) ; % delete insertion
        for n=find(reference.seqvect(5,:)>1.1)
            reference.seqvect = [ reference.seqvect(:,1:n) [x(randperm(4));1] reference.seqvect(:,(n+1):end)] ; % insert deletion
            reference.seqvect(5,n) = 1 ; % reset persistence counter
        end
    end
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    [~,Ali] = swalign(reference.Sequence,true_seq(rndset).Sequence);
    d = seqpdist([Ali(1,:);Ali(3,:)],'Method','p-distance') ;
    fprintf( pid,'%.4f, %d, %d, %d\n', 1-d, L, length(reference.seqvect), length(true_seq(rndset).Sequence)); 
end
fclose(pid);
