function [] = plot_dtwnucleo(tree, weight, true_seq, reference, viewalign, fileout)
% Visualization function to generate nucleotide sequences and plot phylogenetic tree
% tree is a evolving tree structure output by ETDTWrec function

    if nargin<6
        fileout=0 ;
        if nargin<5
            viewalign=0 ;
            if nargin<4
                reference=[];
                if nargin<3
                    true_seq=[];
                end
            end
        end
    end

    % in case reference has extra fields
    if isstruct(reference)
        ref.Header = reference.Header ;
        ref.Sequence = reference.Sequence ;
        ref.seqvect = reference.seqvect(1:4,:) ; % get rid of persistence vector
    else
        ref=[];
    end
    
    % if a tree is input then only consider the tips
    if (numel(tree)>1)
        tips = setdiff(1:length(tree),tree);
    else % otherwise consider all weight matrices
        tips = 1:numel(weight) ;
    end
    
    if iscell(weight)
        weight_seq=[];
        s=1;
        for n=tips
            weight_seq(s).Header = ['weight' num2str(n)] ;
            weight_seq(s).Sequence = mat2nucleo(weight{n}) ;
            weight_seq(s).seqvect = weight{n}(1:4,:) ; % get rid of persistence vector
            s=s+1;
        end
    else
        weight_seq=[];
    end
    
    % make sure the data structure dimension are consistent
    if size(weight_seq,2)>size(weight_seq,1), weight_seq=weight_seq'; end
    if size(true_seq,2)>size(true_seq,1), true_seq=true_seq'; end
    
    if (viewalign>0) % multiple seq alignment
        seqalignviewer(multialign([true_seq ; ref ; weight_seq])) ;
    end
    
    phytree = seqlinkage(seqpdist([true_seq ; ref ; weight_seq],'Indels','pairwise-delete'),'average',[true_seq ; ref ; weight_seq]); 
    phytreeviewer(phytree);
    
    if fileout~=0
        fastawrite(fileout,[true_seq ; ref ; weight_seq]);
    end
    
end