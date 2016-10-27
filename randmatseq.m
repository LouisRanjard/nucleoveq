function [ randseq ] = randmatseq(len)
% return a random sequence of nucleotides of length len,
% encoded as 4 element vectors
% each base is at same frequency (roughly if mod(len,4)>0)

    a = diag(ones(1,4)) ;
    b = repmat(a,1,len) ;
    randseq = b(:,randperm(len)) ;

end
