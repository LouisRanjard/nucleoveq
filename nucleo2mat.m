function [ matseq ] = nucleo2mat(nucleoseq, persistence)
% convert nucleotide sequence to vectors, if persistence

if nargin<2
    persistence=0;
end

if persistence==0
    matseq = zeros(4,length(nucleoseq)) ;
    for nt=1:length(nucleoseq)
        switch nucleoseq(nt) ;
            case 'N'
                matseq(:,nt) = [ 1 1 1 1 ];
            case 'A'
                matseq(:,nt) = [ 1 0 0 0 ];
            case 'C'
                matseq(:,nt) = [ 0 1 0 0 ];
            case 'T'
                matseq(:,nt) = [ 0 0 1 0 ];
            case 'G'
                matseq(:,nt) = [ 0 0 0 1 ];
        end
    end
else
    matseq = zeros(5,length(nucleoseq)) ;
    for nt=1:length(nucleoseq)
        switch nucleoseq(nt) ;
            case 'N'
                matseq(:,nt) = [ 1 1 1 1 1 ];
            case 'A'
                matseq(:,nt) = [ 1 0 0 0 1 ];
            case 'C'
                matseq(:,nt) = [ 0 1 0 0 1 ];
            case 'T'
                matseq(:,nt) = [ 0 0 1 0 1 ];
            case 'G'
                matseq(:,nt) = [ 0 0 0 1 1 ];
        end
    end
end
    
end
