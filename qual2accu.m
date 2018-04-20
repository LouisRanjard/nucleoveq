function [ accu, phred] = qual2accu(qual,offset)
%
% convert fastq encoded quality to Phred quality scores to accuracy
%
% offset is:
% 33 for Sanger or Illumina 1.8+
% 64 for Solexa+64, Illumina 1.3+, Illumina 1.5+ 
%

if nargin<2
    offset=33;
end

phred = double(qual)-offset ;
accu = 10.^(-phred/10) ;
