function [] = plot_dtwnucleo(tree, weight, true_seq, reference, viewalign, fileout)
% Visualization function to generate nucleotide sequences and plot phylogenetic tree

    if nargin<6
        fileout=0 ;
        if nargin<5
            viewalign=0 ;
        end
    end

    % in case reference has extra fields
    ref.Header = reference.Header ;
    ref.Sequence = reference.Sequence ;
    ref.seqvect = reference.seqvect ;
    
    tips = setdiff(1:length(tree),tree);
    weight_seq=[];
    s=1;
    for n=tips
        weight_seq(s).Header = ['weight' num2str(n)] ;
        weight_seq(s).Sequence = mat2nucleo(weight{n}) ;
        weight_seq(s).seqvect = weight{n} ;
        s=s+1;
    end
    
    if (viewalign>0) % multiple seq alignment
        seqalignviewer(multialign([true_seq ; ref ; weight_seq'])) ;
    end
    
    phytree = seqlinkage(seqpdist([true_seq ; ref ; weight_seq']),'average',[true_seq ; ref ; weight_seq']); 
    phytreeviewer(phytree);
    
    if fileout~=0
        fastawrite(fileout,[true_seq ; ref ; weight_seq']);
    end
    
end