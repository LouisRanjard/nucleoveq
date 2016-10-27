function [ ] = readcoverage(reads)
% get the coverage for a set of reads whose positions are included in Header
%
% reads = short read sequences

    % read length
    readlen = abs(mean(arrayfun(@(x) length(x.seqvect), reads))) ;

    % maximum haplotype length
    len = 0 ;
    for n=1:numel(reads)
        position=str2double(reads(n).Header(strfind(reads(n).Header,'pos')+3:end)) ;
        if position>len
            len=position;
        end
    end
    
    % calculate coverage
    coverage=zeros(1,len);
    for n=1:numel(reads)
        position=str2double(reads(n).Header(strfind(reads(n).Header,'pos')+3:end)) ; 
        for m=position:(position+readlen-1)
            coverage(m)=coverage(m)+1;
        end
    end
    
    % plot
    bar(coverage);

end
