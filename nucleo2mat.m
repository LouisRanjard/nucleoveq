function [ matseq ] = nucleo2mat(nucleoseq, qualiti, persistence)
% convert nucleotide sequence to vectors, if persistence add a vector of 1s

if nargin<3
    persistence=0;
    if nargin<2
        qualiti=zeros(1,length(nucleoseq));
    end
end

if (length(qualiti)~=length(nucleoseq))
    error('nucleo2mat(): Quality and Sequence lengths mismatch');
end

matseq = zeros(4,length(nucleoseq)) ; % [ A C T G ]
for nt=1:length(nucleoseq)
    switch nucleoseq(nt)
        case {'N','-','.'}
            matseq(:,nt) = [ .25 .25 .25 .25 ];
        case 'A'
            matseq(:,nt) = [ 1-qualiti(nt) qualiti(nt)/3 qualiti(nt)/3 qualiti(nt)/3 ];
        case 'C'
            matseq(:,nt) = [ qualiti(nt)/3 1-qualiti(nt) qualiti(nt)/3 qualiti(nt)/3 ];
        case {'T','U'}
            matseq(:,nt) = [ qualiti(nt)/3 qualiti(nt)/3 1-qualiti(nt) qualiti(nt)/3 ];
        case 'G'
            matseq(:,nt) = [ qualiti(nt)/3 qualiti(nt)/3 qualiti(nt)/3 1-qualiti(nt) ];
        case 'R'
            matseq(:,nt) = [ (1-qualiti(nt))/2 qualiti(nt)/2 qualiti(nt)/2 (1-qualiti(nt))/2 ]; % A,G
        case 'Y'
            matseq(:,nt) = [ qualiti(nt)/2 (1-qualiti(nt))/2 (1-qualiti(nt))/2 qualiti(nt)/2 ]; % C,T
        case 'S'
            matseq(:,nt) = [ qualiti(nt)/2 (1-qualiti(nt))/2 qualiti(nt)/2 (1-qualiti(nt))/2 ]; % G,C
        case 'W'
            matseq(:,nt) = [ (1-qualiti(nt))/2 qualiti(nt)/2 (1-qualiti(nt))/2 qualiti(nt)/2 ]; % A,T
        case 'K'
            matseq(:,nt) = [ qualiti(nt)/2 qualiti(nt)/2 (1-qualiti(nt))/2 (1-qualiti(nt))/2 ]; % G,T
        case 'M'
            matseq(:,nt) = [ (1-qualiti(nt))/2 (1-qualiti(nt))/2 qualiti(nt)/2 qualiti(nt)/2 ]; % A,C
        case 'B'
            matseq(:,nt) = [ qualiti(nt) (1-qualiti(nt))/3 (1-qualiti(nt))/3 (1-qualiti(nt))/3 ]; % C,G,T
        case 'D'
            matseq(:,nt) = [ (1-qualiti(nt))/3 qualiti(nt) (1-qualiti(nt))/3 (1-qualiti(nt))/3 ]; % A,G,T
        case 'H'
            matseq(:,nt) = [ (1-qualiti(nt))/3 (1-qualiti(nt))/3 (1-qualiti(nt))/3 qualiti(nt) ]; % A,C,T
        case 'V'
            matseq(:,nt) = [ (1-qualiti(nt))/3 (1-qualiti(nt))/3 qualiti(nt) (1-qualiti(nt))/3 ]; % A,C,G
    end
end

if persistence~=0
    matseq = [matseq; ones(1,length(nucleoseq))] ;
end
    
end
