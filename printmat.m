function [ ] = printmat(matseq,a,b)
% Print a sequence of vector converting to nucleotide sequence

    if nargin<2
        a=1;
        b=length(matseq);
    end

%     tmp=matseq(:,a:b);
%     T=array2table(tmp);
%     P=strtrim(cellstr(num2str([a:b]'))');
%     T.Properties.VariableNames=strcat( num2cell(mat2nucleo(tmp)), P) ;
%     disp(T);
    
    L=length(num2str(b)); % get the number of digits required
    for n=a:b
        fprintf(1,['%0' num2str(L) 'd\t'],n);
    end
    fprintf(1,'\n');
    for n=a:b
        fprintf(1,'<strong>%s</strong>\t',mat2nucleo(matseq(1:4,n)));
    end
    fprintf(1,'\n');
    
    for m=1:length(matseq(:,1))
        for n=a:b
            fprintf(1,'%.4f\t',matseq(m,n));
        end
        fprintf(1,'\n');
    end

end