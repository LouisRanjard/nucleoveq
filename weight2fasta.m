function [ weight_seq ] = weight2fasta(treeH,weightH,writef)
% generate fasta sequences from weight sequences
% write them to fasta file if writef==1
%
% treeH = evolving classification tree (ETDTWrec) 
% weightH = weight matrices
%
    if (numel(treeH)>0)
        tips = setdiff(1:length(treeH),treeH) ;
    else
        tips = 1:numel(weightH) ;
    end
    weight_seq=[];
    s=1;
    for n=tips
        weight_seq(s).Header = ['weight' num2str(n)] ;
        weight_seq(s).seqvect = weightH{n} ;
        weight_seq(s).Sequence = mat2nucleo(weightH{n}) ;
        s=s+1;
    end
    
    % write to file if needed
    if nargin>2 && numel(writef)>0
        fastawrite(writef,weight_seq);
    end
    
end