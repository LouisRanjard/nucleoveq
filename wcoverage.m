function [wcoverage, varcov, shannonEnt ] = wcoverage(weight,readlen,rposition,plotting)
% plot the read coverage and the entropy for a input weight matrix
%
% weight: 4 columns sequence of nucleotides
% readlen: read length
% rposition: alignment position of a set of reads to the weight sequence
%
% wcoverage: coverage of each position in the sequence (=number of reads covering it)
% varcov: normalised coverage variance
% shannonEnt: Shannon entropy of the 4 columns sequence
%
% ** Warning: coverage is calculated by ignoring the indels in the alignment of reads to reference **

    wcoverage=nan ; 
    varcov=nan ; 
    shannonEnt=nan ;
    
    if numel(rposition)==0
        return;
    end

    if nargin<4
        plotting=false;
    end
    
    % coverage
    wcoverage=zeros(1,size(weight,2));
    for n=1:length(rposition)
        for m=max([ 1 (rposition(n)-readlen+1) ]):rposition(n)
            wcoverage(m)=wcoverage(m)+1;
        end
    end
    
    % normalised coverage variance
    varcov = var(wcoverage/norm(wcoverage)) ;
    
    % plot with entropy
    if plotting
        [ ~, shannonEnt ] = shannonEntropy(weight) ;
        [hAx,hLine1,hLine2] = plotyy(1:length(shannonEnt), wcoverage, 1:length(shannonEnt), shannonEnt, @bar, @line); % R2015b
        %xlabel('Position');
        %ylabel(hAx(1),'Coverage'); % left y-axis
        %ylabel(hAx(2),'Entropy');  % right y-axis
        set(hAx,'xtick',[]); set(hAx,'xticklabel',[]);
        set(hAx(1),'ytick',[]); set(hAx(1),'yticklabel',[]);
        set(hAx(2),'ytick',[]); set(hAx(2),'yticklabel',[]);
        set(hLine2,'LineWidth',2);
        set(hLine1,'EdgeColor',[.7 .7 .7]);
        set(hLine1,'FaceColor',[.7 .7 .7]);
        set(hAx,'XLim',[0 length(shannonEnt)]) ;
    %     yyaxis left ; R2016a
    %     bar(coverage) ;
    %     hold on ;
    %     yyaxis right ;
    %     plot(shannonEnt) ;
    else
        shannonEnt=nan ;
    end

end