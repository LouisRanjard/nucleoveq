function [ mutmatseq ] = mutatematseq(matseq,percent,percentindels)
% return a mutated sequence inserting "percent" substitutions
% and "percentindels" insertions and deletions 
% 
% ex: mutatematseq(matseq,0.1) for 90% similarity 

    if sum( reshape(matseq(1:4,:),[],1)>0 & reshape(matseq(1:4,:),[],1)<1 )>0
        error('need fully defined matrix (0 and 1 only)');
    end

    if nargin<3, percentindels=0; end

    mutmatseq = matseq ;

    % SUBSTITUTIONS
    subpos = datasample(1:length(matseq),round(length(matseq)*percent),'Replace',false) ;
    for n=1:numel(subpos) % each position to be mutated
        zero = find(matseq(:,subpos(n))==0) ; % row positions that ==0 and can be changed
        mutmatseq(1:4,subpos(n)) = zeros(4,1) ; % set the position to zero
        m = randi([1 3],1,1) ;
        mutmatseq(zero(m),subpos(n)) = 1 ;
        %if sum([ matseq(:,subpos(n)) - mutmatseq(:,subpos(n)) ]==0)==0
        %    disp(subpos(n));
        %end
    end
    
    % INDELS
    indelpos = datasample(1:(length(matseq)-1),round((length(matseq)-1)*percentindels),'Replace',false) ;
    for n=1:numel(indelpos) % each position to be mutated
        if random('unif',0,1)<0.5 % DELETION
            mutmatseq = [mutmatseq(:,1:(indelpos(n)-1)) mutmatseq(:,(indelpos(n)+1):end)] ;
        else % INSERTION
            newpos = zeros(4,1) ;
            newpos(randi([1 4],1,1)) = 1 ; % choose nucleotide to insert at random
            if size(mutmatseq,1)>4
                newpos=[newpos; 1]; % add persistence vector
            end
            mutmatseq = [mutmatseq(:,1:indelpos(n)) newpos mutmatseq(:,(indelpos(n)+1):end)] ;
        end
    end
    
    if sum(sum(mutmatseq(1:4,:),1)~=1)>0 % check if each column (locus) sums to 1
        error('sequence contains undefined positions');
    end
    
end
