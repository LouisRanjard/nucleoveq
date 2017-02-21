function [ matseqRC ] = reversecomplement(matseq)
% return reverse complement of matseq sequence
% A -> T
% C -> G
% T -> A
% G -> C

    matseqRC = matseq([3 4 1 2],end:-1:1) ;

end