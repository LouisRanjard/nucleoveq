function [ mutmatseq ] = mutatematseq(matseq,percent)
% return a mutated sequence inserting "percent" mutations
% ex: mutatematseq(matseq,0.1) for 90% similarity 

    mutmatseq = matseq ;
    mutationpos = datasample(1:length(matseq),length(matseq)*percent,'Replace',false) ;

    for n=1:numel(mutationpos) % each position to be mutated
        zero = find(matseq(:,mutationpos(n))==0) ; % row positions that ==0 and can be changed
        mutmatseq(:,mutationpos(n)) = zeros(4,1) ; % set the position to zero
        m = randi([1 3],1,1) ;
        mutmatseq(zero(m),mutationpos(n)) = 1 ;
        %if sum([ matseq(:,mutationpos(n)) - mutmatseq(:,mutationpos(n)) ]==0)==0
        %    disp(mutationpos(n));
        %end
    end
    
    if sum(sum(mutmatseq,1)~=1)>0 % check if each column (locus) sums to 1
        error('sequence contains undefined positions');
    end
end
