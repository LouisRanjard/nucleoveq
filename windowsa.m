function [ wcovera, wentrop, xtl ] = windowsa(matseq,siz,ovlap,BMU,plotting)
% compute statistics of input matrix matseq using rectangular window
% analysis of size "siz" overlapping with "ovlap"*100 percent
%
% return:
% wcovera: coverage at each window as the number of reads which mapping starts inside
% wentrop: entropy of each window
% xtl: XTickLabel with starting position of each window in the original matseq

    wcovera = zeros() ;
    wentrop = zeros() ;
    
    n = 1 ;
    lag = round(siz*(1-ovlap)) ;
    if lag==0, lag=1; end % at minimum, move every position in the alignment
    xtl = 1:lag:size(matseq,2) ;
    for wstart=xtl
       wend = min( [wstart+siz ; size(matseq,2)] ) ;
       %wcovera(n) = sum( BMU(:,3)>=wstart & BMU(:,3)<=wend ) ; % all the reads which mapping starts somewhere within the window
       wcovera(n) = sum( BMU(:,3)<=wstart & BMU(:,4)>=wend ) ; % all the reads which mapping starts somewhere within the window
       if nargout>1
         wentrop(n) = shannonEntropy(matseq(1:4,wstart:wend)) ;
       end
       n = n + 1 ;
    end
    
    if plotting==1
        figure;
        yyaxis left
        bar(1:length(wcovera),wcovera) ;
        yyaxis right
        plot(1:length(wentrop),wentrop) ;
    end
    
end

