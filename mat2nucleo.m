function [ nucleoseq ] = mat2nucleo(matseq)
% convert nucleotide vector sequence to fasta

    nucleoseq = repmat(char(0),1,length(matseq));
    for nt=1:length(matseq)
        %[~, nucleo] = max(matseq(1:4,nt)) ;
        nucleo = find(matseq(1:4,nt)==max(matseq(1:4,nt))) ;
        if numel(nucleo)>1
            nucleoseq(nt) = 'N' ;
        else
            switch nucleo ;
                case 1
                    nucleoseq(nt) = 'A' ;
                case 2
                    nucleoseq(nt) = 'C' ;
                case 3
                    nucleoseq(nt) = 'T' ;
                case 4
                    nucleoseq(nt) = 'G' ;
            end
        end
    end
    
    %{
    for ds=1:numel(dnaseq)
    dnaseq(ds).Sequence=[];
    for nt=1:length(dnaseq(ds).seqvect)
        if ( numel(find(dnaseq(ds).seqvect(nt,1:4)==1)) == 1 )
            switch find(dnaseq(ds).seqvect(nt,1:4)==1) ;
                case 1
                    dnaseq(ds).Sequence = [ dnaseq(ds).Sequence; 'A' ];
                case 2
                    dnaseq(ds).Sequence = [ dnaseq(ds).Sequence; 'C' ];
                case 3
                    dnaseq(ds).Sequence = [ dnaseq(ds).Sequence; 'T' ];
                case 4
                    dnaseq(ds).Sequence = [ dnaseq(ds).Sequence; 'G' ];
            end
        elseif ( numel(find(dnaseq(ds).seqvect(nt,1:4)==1)) == 2 )
            if ( dnaseq(ds).seqvect(nt,1)==1 && dnaseq(ds).seqvect(nt,4)==1 )
                dnaseq(ds).Sequence = [ dnaseq(ds).Sequence; 'R' ];
            elseif ( dnaseq(ds).seqvect(nt,2)==1 && dnaseq(ds).seqvect(nt,3)==1 )  
                dnaseq(ds).Sequence = [ dnaseq(ds).Sequence; 'Y' ]; 
            else
                dnaseq(ds).Sequence = [ dnaseq(ds).Sequence; 'N' ];
            end
        else
            dnaseq(ds).Sequence = [ dnaseq(ds).Sequence; 'N' ];
        end
    end
    end
    %}

end

