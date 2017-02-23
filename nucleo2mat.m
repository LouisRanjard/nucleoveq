function [ matseq ] = nucleo2mat(nucleoseq, qualiti, persistence)
% convert nucleotide sequence to vectors, if persistence

if nargin<3
    persistence=0;
    if nargin<2
        qualiti=zeros(1,length(nucleoseq));
    end
end

if (length(qualiti)~=length(nucleoseq))
    error('nucleo2mat(): Quality and Sequence lengths mismatch');
end

if persistence==0
    matseq = zeros(4,length(nucleoseq)) ;
    for nt=1:length(nucleoseq)
        switch nucleoseq(nt)
            case {'N','-'}
                matseq(:,nt) = [ .25 .25 .25 .25 ];
            case 'A'
                matseq(:,nt) = [ 1-qualiti(nt) qualiti(nt)/3 qualiti(nt)/3 qualiti(nt)/3 ];
            case 'C'
                matseq(:,nt) = [ qualiti(nt)/3 1-qualiti(nt) qualiti(nt)/3 qualiti(nt)/3 ];
            case 'T'
                matseq(:,nt) = [ qualiti(nt)/3 qualiti(nt)/3 1-qualiti(nt) qualiti(nt)/3 ];
            case 'G'
                matseq(:,nt) = [ qualiti(nt)/3 qualiti(nt)/3 qualiti(nt)/3 1-qualiti(nt) ];
        end
    end
else
    matseq = zeros(5,length(nucleoseq)) ;
    for nt=1:length(nucleoseq)
        switch nucleoseq(nt)
            case {'N','-'}
                matseq(:,nt) = [ .25 .25 .25 .25 1 ];
            case 'A'
                matseq(:,nt) = [ 1-qualiti(nt) qualiti(nt)/3 qualiti(nt)/3 qualiti(nt)/3 1 ];
            case 'C'
                matseq(:,nt) = [ qualiti(nt)/3 1-qualiti(nt) qualiti(nt)/3 qualiti(nt)/3 1 ];
            case 'T'
                matseq(:,nt) = [ qualiti(nt)/3 qualiti(nt)/3 1-qualiti(nt) qualiti(nt)/3 1 ];
            case 'G'
                matseq(:,nt) = [ qualiti(nt)/3 qualiti(nt)/3 qualiti(nt)/3 1-qualiti(nt) 1 ];
        end
    end
end
    
end
