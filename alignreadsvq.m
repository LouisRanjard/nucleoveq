function [weight_seq,rpositionr,rdist] = alignreadsvq(reads1,reads2,reffasta,options)
%% wrapper function to align reads to a reference using vector quantization
%
% options.setsize : number of weight matrices to recosntruct
% options.minproportion : minimum proportion of weight matrices (if not Uniform) to use to compute learning rates
%
    % initialise seed at random
    fprintf(1,'%s - Initialising...\n', datestr(now)) ;
    rng('shuffle');
    
    pe=1;
    if isempty(reads2)
        pe=0;
    end
    
    % how long to consider before and after mapping position when MAPPING file is provided
    try width=options.width; catch width=+Inf ; end
    
    % indels parameters
    try indel_cost=options.indel_cost; catch indel_cost=1 ; end
    try indel_weight=options.indel_weight; catch indel_weight=2 ; end
    
    % encode the reads
    if isfield(reads1,'Quality')
        for m=1:numel(reads1)
            reads1(m).seqvect=nucleo2mat(reads1(m).Sequence,qual2accu(reads1(m).Quality)) ;
            if pe, reads2(m).seqvect=nucleo2mat(reads2(m).Sequence,qual2accu(reads2(m).Quality)) ;end
        end
    else
        for m=1:numel(reads1)
            reads1(m).seqvect=nucleo2mat(reads1(m).Sequence) ;
            if pe, reads2(m).seqvect=nucleo2mat(reads2(m).Sequence) ;end
        end
    end
    NREAD=numel(reads1) ;
    % for each read1 we need a matching read2, control tbd here

    % load reference
    if nargin<3
      reffasta='/home/louis/Documents/Projects/Pooling3/Macropodidae/EasternGrey_NC027424_Amplicons.fasta';
    end
    %reference = ref_fasta(1) ; % sample 10 is Amplicon 1
    try reference=options.reference; catch reference=fastaread(reffasta) ; end
    if ~isfield(reference,'seqvect')
       reference.seqvect = nucleo2mat(reference.Sequence) ;
    end
    if size(reference.seqvect,1)<6
       reference.seqvect = [reference.seqvect ; ...
                          ones(1,length(reference.Sequence)) ; ... % add persistence vector
                          zeros(1,length(reference.Sequence))] ; % add coverage vector
    end

    % Number of weight matrices to project the reads on
    try numw=options.setsize; catch numw=3; end
    % mapping of the reads to the reference
    try mapfile=options.mapping; catch mapfile=[]; end
    if numel(mapfile)>0 && numw>0
        warning('SAM file positioning is not implemented yet for multiple weights -> SAM ignored');
        mapfile=[];
    end
    % trim the ends of the reconstructed sequence when the coverage is low (set to 0 to not do any trimming)
    try trimend=options.trimend; catch trimend=0.05; end
    try w1=options.weight; catch w1=0; end
    if (w1>0)
        cw=0; % use the input weight value
    else
        cw=1; % compute the weight (learning rate) depending on distance
    end
    
    % NOT USED: root weight initialisation: align all reads once to the reference to create initial weight
    N = numel(reads1) + numel(reads2) ; % total number of reads
    if pe % read lengths
        rlen=cellfun(@(x) length(x), {reads1(:).Sequence reads2(:).Sequence});
    else
        rlen=cellfun(@(x) length(x), {reads1(:).Sequence});
    end
    Lr=sum(rlen)/length(rlen); % average Length of reads
    LR=length(reference.Sequence); % Length of Reference
    cover1 = (N*Lr)/LR ; % expected coverage of the reference matrix
    %cover1 = 15; % needs to be estimated prior?
    %w1 = nthroot(0.001,(cover1/2)-1) ; % doesn't work when covrage ~2
    %w1=0.5;
    %fprintf( 1,'cover1=%.2f, w1=%.5f\n',cover1,w1 );
    fprintf( 1,'coverage=%.2f\n',cover1 );
    try minP=options.minproportion; catch minP=1/max(1,numw); end % by default average coverage of each weight matrix
    cover2 = cover1*minP ; % use different options.minproportion to e.g. use the lowest coverage for Pooling3 experiments ( 12.5% )
    %cover2 = cover1*0.0625 ; % using half the lowest coverage for Pooling3 experiments ( 12.5% / 2 )
    %cover2 = cover1*0.01 ; % speed up learning?
    %w2 = nthroot(0.001,(cover2/2)-1) ; % doesn't work when covrage ~2
    w2 = w1 ;
    %fprintf( 1,'cover2=%.2f, w2=%.5f\n',cover2,w2 );
    
    if numel(mapfile)~=0
        fprintf(1,'%s - Retrieving read mapping positions...\n', datestr(now)) ;
        mapping = samread(mapfile) ;
        cells = struct2cell(mapping) ; % converts struct to cell matrix
        % get the mapping position of each read
        mapObj = containers.Map(cells(1,:),cells(4,:)) ; % build a "hash table"
        L_all = length(reads1)+length(reads2) ;
        reads1m = reads1(isKey(mapObj,{reads1.Header})) ; % only keep mapped reads
        if pe
            reads2m = reads2(isKey(mapObj,{reads2.Header})) ;
        end
        reads1u = reads1(~isKey(mapObj,{reads1.Header})) ; % unmapped reads
        if pe, reads2u = reads2(~isKey(mapObj,{reads2.Header})) ;end
        p1 = values(mapObj,{reads1m.Header}) ; % retrieve mapping position
        r1b = squeeze(struct2cell(reads1m));
        r1b = [r1b;p1];
        reads1m = cell2struct(r1b,[fieldnames(reads1m)' 'Position' ],1);
        [reads1u(:).Position] = deal(NaN) ; % set to '0' the position of unmapped reads
        if pe
            p2 = values(mapObj,{reads2m.Header}) ; % retrieve mapping position
            r2b = squeeze(struct2cell(reads2m));
            r2b = [r2b;p2];
            reads2m = cell2struct(r2b,[fieldnames(reads2m)' 'Position' ],1);
            [reads2u(:).Position] = deal(NaN) ;
        end
        if pe
            %L_mapped = length(reads1m)+length(reads2m) ;
            L_mapped = sum([reads1m.Position]>0) + sum([reads2m.Position]>0) ;
        else
            %L_mapped = length(reads1m) ;
            L_mapped = sum([reads1m.Position]>0) ;
        end
        fprintf(1,'%d reads are mapped out of %d reads\n',L_mapped,L_all);
        % ordering the reads by position
        %fprintf(1,'%s - Ordering read mapping positions...\n', datestr(now)) ;
        % order reads by Position in mapping file (would be much faster using awk...)
        %cells = struct2cell(mapping); % converts struct to cell matrix
        %sortvals = cells(4,:);        % gets the values of just the Position
        %mat = cell2mat(sortvals);     % converts values to a matrix
        %mat = squeeze(mat);           % removes the empty dimensions for a single vector
        %[~,ix] = sort(mat);           % sorts the vector of values
        %mapping = mapping(ix);        % rearranges the original array
        % do same order for reads
        %cells = struct2cell(mapping); % converts position ordered struct to cell matrix
        %names = cells(1,:);           % get mapping read names only
        %x = cellfun(@(x) strfind(x,' 1:'), names, 'UniformOutput', false); % find reads 1 [SLOW!]
        %y = cellfun(@numel,x);
        %z = names(logical(y));        % get the ordered names of reads1 only
        %rindex = cellfun(@(x) find(strcmp({reads1.Header},x)), names);
        % join the two structures and ignore the pairing (simpler because
        % some pairs may have been broken when only one of the two can be
        % mapped
        if pe
            reads = [reads1m; reads1u'; reads2m; reads2u'];
        else
            reads = [reads1m; reads1u'];
        end
        %rindex = randperm(numel(reads)) ; % random order
        for i=1:numel(reads) % set to NaN the reads that were unmapped and cast Position to DOUBLE (from UNIT32) because only DOUBLE can deal with NaN
            if reads(i).Position==0
                reads(i).Position = NaN ;
            else
                reads(i).Position = double(reads(i).Position) ;
            end
        end
        [~,rindex] = sort([reads(:).Position]); % sort by mapping Position (automatically ignore NaN and put them at the end of sorted index list)
        rpositionr = zeros(numel(reads),1) ; % position of the alignment of the first read in each pair
        rdist = zeros(numel(reads),1);
        fprintf(1,'%s - Building root matrix...\n', datestr(now)) ;
        counter=1;
        for r=rindex
%fprintf(1,'%d %d\n',size(reference.seqvect));
%plotnucleoveq(reference) ; title([num2str(counter) '/' num2str(numel(reads))]); drawnow; counter=counter+1;
            % get mapping position
            if isnan(reads(r).Position)
                a = 1;
                b = size(reference.seqvect,2);
            else
                a = max( [ 1 reads(r).Position-width ] ) ;
                b = min( [ size(reference.seqvect,2) a+length(reads(r).Sequence)+width ] ) ;
            end
            % align read 1
            [ distF, seqvectF, align_endF ] = DTWaverage( reference.seqvect(:,a:b), reads(r).seqvect, 1, w1, 0, 1, cw, indel_cost, indel_weight ) ;
%fprintf(1,'%f  %d\n',distF,sum(seqvectF(5,:)~=1)) ;
            RC = reversecomplement(reads(r).seqvect);
            [ distR, seqvectR, align_endR ] = DTWaverage( reference.seqvect(:,a:b), RC, 1, w1, 0, 1, cw, indel_cost, indel_weight ) ;
%fprintf(1,'%f  %d\n',distR,sum(seqvectR(5,:)~=1)) ;
            if distF<distR
              reference.seqvect = [ reference.seqvect(:,1:a-1) seqvectF reference.seqvect(:,b+1:end) ] ;
              rpositionr(r) = align_endF + a ;
              rdist(r) = distF ;
            else
              reads(r).seqvect = RC ; % only keep FORWARD sequence because there is no direction check in ETDTWrec
              if isfield(reads,'Quality'), reads(r).Quality = reads(r).Quality(end:-1:1) ; end
              reads(r).Sequence = mat2nucleo(reads(r).seqvect) ;
              reference.seqvect = [ reference.seqvect(:,1:a-1) seqvectR reference.seqvect(:,b+1:end) ] ;
              rpositionr(r) = align_endR + a ;
              rdist(r) = distR ;
            end
        end
    else % no mapping file provided
        % the reads are paired, align one after the other
        try rindex=options.rindex; catch rindex = randperm(numel(reads1)) ; end
        rpositionr = zeros(numel(reads1),1) ; % position of the alignment of the first read in each pair
        if pe
            rdist = zeros(numel(reads1)+numel(reads2),1);
        else
            rdist = zeros(numel(reads1),1);
        end
        fprintf(1,'%s - Building root matrix...\n', datestr(now)) ;
        a = 1 ;
        c = 1 ;
        for r=rindex
            %fprintf(1,'%d %d\n',size(reference.seqvect));
            % align read 1
            b = size(reference.seqvect,2) ;
            [ distF, seqvectF, align_endF ] = DTWaverage( reference.seqvect(:,a:b), reads1(r).seqvect, 1, w1, 0, 1, cw, indel_cost, indel_weight ) ;
            RC = reversecomplement(reads1(r).seqvect);
            [ distR, seqvectR, align_endR ] = DTWaverage( reference.seqvect(:,a:b), RC, 1, w1, 0, 1, cw, indel_cost, indel_weight ) ;
            if distF<distR
              reference.seqvect = [ reference.seqvect(:,1:a-1) seqvectF reference.seqvect(:,b+1:end) ] ;
              rpositionr(r) = align_endF + a ;
              rdist(r) = distF ;
            else
              reads1(r).seqvect = RC ; % only keep FORWARD sequence because there is no direction check in ETDTWrec
              if isfield(reads1,'Quality'), reads1(r).Quality = reads1(r).Quality(end:-1:1) ; end
              reads1(r).Sequence = mat2nucleo(reads1(r).seqvect) ;
              reference.seqvect = [ reference.seqvect(:,1:a-1) seqvectR reference.seqvect(:,b+1:end) ] ;
              rpositionr(r) = align_endR + a ;
              rdist(r) = distR ;
            end
            if pe % align read 2
                d = size(reference.seqvect,2) ;
                [ distF, seqvectF, align_endF ] = DTWaverage( reference.seqvect(:,c:d), reads2(r).seqvect, 1, w1, 0, 1, cw, indel_cost, indel_weight ) ;
                RC = reversecomplement(reads2(r).seqvect);
                [ distR, seqvectR, align_endR ] = DTWaverage( reference.seqvect(:,c:d), RC, 1, w1, 0, 1, cw, indel_cost, indel_weight ) ;
                if distF<distR
                  reference.seqvect = [ reference.seqvect(:,1:c-1) seqvectF reference.seqvect(:,d+1:end) ] ;
                  rdist(r+NREAD) = distF ;
                  if align_endF<rpositionr(r) % the read 2 in pair aligns first in the sequence (before read1) -> swap them
                      readtmp=reads1(r); reads1(r)=reads2(r); reads2(r)=readtmp;
                  end 
                else
                  reads2(r).seqvect = RC ; % only keep FORWARD sequence because there is no direction check in ETDTWrec
                  if isfield(reads2,'Quality'), reads2(r).Quality = reads2(r).Quality(end:-1:1) ; end
                  reads2(r).Sequence = mat2nucleo(reads2(r).seqvect) ;
                  reference.seqvect = [ reference.seqvect(:,1:c-1) seqvectR reference.seqvect(:,d+1:end) ] ;
                  rdist(r+NREAD) = distR ;
                  if align_endR<rpositionr(r) % the read 2 in pair aligns first in the sequence (before read1) -> swap them
                      readtmp=reads1(r); reads1(r)=reads2(r); reads2(r)=readtmp;
                  end
                end
            end
        end
%         % realign using rposition
%         [~,rindex]=sort(rpositionr);
%         for r=rindex'
%           %plotnucleoveq(reference) ; title([num2str(counter) '/' num2str(numel(reads))]); drawnow; counter=counter+1;
%           [ ~, reference.seqvect ] = DTWaverage( reference.seqvect, reads1(r).seqvect, 1, w1, 0, 1 ) ;
%           [ ~, reference.seqvect ] = DTWaverage( reference.seqvect, reads2(r).seqvect, 1, w1, 0, 1 ) ;
%         end
    end
    
    reference.Sequence = mat2nucleo(reference.seqvect) ;
    
    if trimend~=0
        fprintf(1,'%s - Trimming ends...\n', datestr(now)) ;
        %%%%%%%%%%%% trim ends only with no coverage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        deb = find(reference.seqvect(6,:)>0,1) ;
        fin = find(reference.seqvect(6,:)>0,1,'last') ;
        reference.seqvect = reference.seqvect(:,deb:fin) ;
        reference.Sequence = reference.Sequence(deb:fin) ;
        %%%%%%%%%%%% trim ends only below threshold %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         covlimit = quantile(reference.seqvect(6,reference.seqvect(6,:)>=1),trimend) ; % only considers position with actual coverage
%         deb = find(reference.seqvect(6,:)>covlimit,1) ;
%         fin = find(reference.seqvect(6,:)>covlimit,1,'last') ;
%         reference.seqvect = reference.seqvect(:,deb:fin) ;
%         reference.Sequence = reference.Sequence(deb:fin) ;
        %%%%%%%%%%% trim all position below coverage limit %%%%%%%%%%%%%%%%%%%%%%%%
%         reference.seqvect = reference.seqvect(:,reference.seqvect(6,:)>covlimit) ;
%         reference.Sequence = reference.Sequence(reference.seqvect(6,:)>covlimit) ;
        %%%%%%%%%%%% remove position with 0 coverage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         reference.seqvect = reference.seqvect(:,reference.seqvect(6,:)>0) ;
%         reference.Sequence = reference.Sequence(reference.seqvect(6,:)>0) ;
        %%%%%%%%%%%% return the longest sequence with coverage >covlimit %%%%%%%%%%%
        % to allow errors one could use https://stackoverflow.com/questions/32456858/find-longest-sequence-of-non-nan-values-but-allow-for-threshold
%         cov = movimean(reference.seqvect(6,:),Lr); % smooth out the coverage first by doing moving mean around average read length
%         X = cov>0; % use 0 because there has been smoothing done
%         contiglengths = diff(find([1,diff(X),1]));
%         contigstarts = find([1,diff(X)]);
%         [~,idx] = max(contiglengths);
%         deb = contigstarts(idx);
%         fin = deb+contiglengths(idx)-1;
%         reference.seqvect = reference.seqvect(:,deb:fin) ;
%         reference.Sequence = reference.Sequence(deb:fin) ;
    end

    if 0 % remove x% of the pairs that align with the worse scores
        toremove = find(rdist>=quantile(rdist,0.99999)) ;
        % remove the matching pairs as well
        toremove = [ toremove-numel(reads1) toremove toremove+numel(reads1)] ;
        toremove = toremove( toremove>0 & toremove<=numel(reads1) ) ; % only keep indexes of the first in pairs
        tokeep = setdiff(1:numel(reads1),toremove) ;
        rpositionr = rpositionr(tokeep) ;
        rdist = rdist(tokeep) ;
        reads1 = reads1(tokeep) ;
        reads2 = reads2(tokeep) ;
    end

%     % test doing an extra realign here ==> does not help
%     fprintf(1,'%s - Realign...\n', datestr(now)) ;
%     for realign=1:5
%         disp(sum(reference.seqvect(6,:))) ;
%         counter=1;
%         for r=rindex
%           plotnucleoveq(reference) ; title([num2str(counter) '/' num2str(numel(reads))]); drawnow; counter=counter+1;
%           [ ~, reference.seqvect ] = DTWaverage( reference.seqvect, reads(r).seqvect, 1, w1, 0, 1 ) ;
%         end
%         reference.Sequence = mat2nucleo(reference.seqvect) ;
%     end
    
    % implement competitive learning when multiple weights are required
    if numw>=1
        if 1 % realign the reads 1 to the reference without updating, just to get the correct distances (rdist) and positions (rpositionr)
            for r=rindex
              [ dist, ~, align_end ] = DTWaverage( reference.seqvect, reads1(r).seqvect, 1, w1, 0, 1, cw, indel_cost, indel_weight ) ;
              rpositionr(r) = align_end ;
              rdist(r) = dist ;
            end
        end

        % choose the read in order to their position in reference rather than random
        %[~,I] = sort(rpositionr) ;
        % order them by distance to reference, small distance first means most frequent haplotypes first
         [~,I] = sort(rdist(1:numel(reads1))) ; % only sort the first half (first aligning read in pair)
        % order them by distance to reference and then position
        %[~,I] = sortrows([rdist rpositionr],[1,2]) ;
        if 0 % consider consecutive overlapping windows, get all reads aligning in it and sort them by distance to reference: DOES NOT HELP MUCH, TO BE CHECKED?
            winsize = 100 ;
            winstep = round(winsize/2) ;
            I=[] ;
            for j = winsize:winstep:length(reference.Sequence)
                %reads_in = find(rpositionr>=max(1,j-winsize) && rpositionr<=j) ; % reads 1 that align within the current window
                rdist_tmp = rdist(1:numel(reads1)) ;
                rdist_tmp(rpositionr<max(1,j-winsize) | rpositionr>j) = +Inf ; % set to +Infinity the distance to these reads
                [rdist_tmp,i] = sort(rdist_tmp) ; % sort these reads according to their distance to reference
                I = [ I ; i(~isinf(rdist_tmp)) ] ;
            end
        end
        % or/and in lexicographic order?
        % random order
        %I = randperm(numel(reads1))' ;

        % reinitialise variable positions
        %[ ~, shannonEnt ] = shannonEntropy(reference.seqvect) ;
        %reference.seqvect(1:4,shannonEnt>.2) = .25 ; % non_zeros_freq=[.95 .05] ; -sum(non_zeros_freq.*log(non_zeros_freq)) : 0.1985

        % Competitive learning
        fprintf(1,'%s - Learning...\n', datestr(now)) ;
        % Declare data structures
        d_r1 = zeros(1,numw) ; % distance read1 to weight matrix
        d_r2 = zeros(1,numw) ; % distance read2 to weight matrix
        end_r1 = zeros(1,numw) ;
        end_r2 = zeros(1,numw) ;
        start_r1 = zeros(1,numw) ;
        start_r2 = zeros(1,numw) ;
        ref.persistent.copy = [ reference.seqvect(1:4,:); ones(1,length(reference.Sequence)) ] ;
        weight(1:numw) = {ref.persistent.copy} ; % WEIGHT INITIALISATION
        seqvect = cell(1,numw) ; % temporary matrices, used to update weights
        BMU = zeros(2*numel(I),4) ; % weight_id, dist, position_start, position_end
        %verbose=0;
        for r=I' % LEARNING
            %[ d1, seqvect1, end1_r1, start1_r1 ] = DTWaverage( weight{1}, reads1(r).seqvect, 1, w2, 0, 1 ) ;
            %[ d2, seqvect2, end2_r1, start2_r1 ] = DTWaverage( weight{2}, reads1(r).seqvect, 1, w2, 0, 1 ) ;
            %[ d3, seqvect3, end3_r1, start3_r1 ] = DTWaverage( weight{3}, reads1(r).seqvect, 1, w2, 0, 1 ) ;
            for n=1:numw
                [ d_r1(n), seqvect{n}, end_r1(n), start_r1(n) ] = DTWaverage( weight{n}, reads1(r).seqvect, 1, w2, 0, 1, cw, indel_cost, indel_weight ) ;
            end
            %if (d1==d2 && d2==d3)
            if pe & all(d_r1==d_r1(1)) % choose where to align the pair read   
                %fprintf(1,'all weight equals\n'); % if distance to all weights is equal (true at least for the first iteration), 
                %  do not do any update with read1 to avoid doing an arbitrary choice which will have an impact with the update performed with read2
                %only do update according to read2, use read2 to choose weight
                % if distance to read2 are all equal again, then choose one weight randomly to start differentiating the weights against each other
                %[ d1, seqvect1, end1_r2, start1_r2 ] = DTWaverage( seqvect1, reads2(r).seqvect, 1, w2, 0, 1 ) ;
                %[ d2, seqvect2, end2_r2, start2_r2 ] = DTWaverage( seqvect2, reads2(r).seqvect, 1, w2, 0, 1 ) ;
                %[ d3, seqvect3, end3_r2, start3_r2 ] = DTWaverage( seqvect3, reads2(r).seqvect, 1, w2, 0, 1 ) ;
                for n=1:numw
                    [ d_r2(n), seqvect{n}, end_r2(n), start_r2(n) ] = DTWaverage( seqvect{n}, reads2(r).seqvect, 1, w2, 0, 1, cw, indel_cost, indel_weight ) ;
                end
                %if (d1==d2 && d2==d3) % all equal, choose at random which weight to update
                if all(d_r2==d_r2(1)) 
%                     x=rand(1);
%                     if x<0.333
%                         weight{1} = seqvect1 ;
%                         BMU(r,:) = [1 d1 start1_r1 end1_r1] ;
%                         BMU(r+NREAD,:) = [1 d1 start1_r2 end1_r2] ;
%                     elseif x<0.666
%                         weight{2} = seqvect2 ;
%                         BMU(r,:) = [1 d2 start2_r1 end2_r1] ;
%                         BMU(r+NREAD,:) = [2 d2 start2_r2 end2_r2] ;
%                     else
%                         weight{3} = seqvect3 ;
%                         BMU(r,:) = [1 d3 start3_r1 end3_r1] ;
%                         BMU(r+NREAD,:) = [3 d3 start3_r2 end3_r2] ;
%                     end
                    x=randsample(numw,1);
                    weight{x} = seqvect{x} ;
                    BMU(r,:) = [x d_r1(x) start_r1(x) end_r1(x)] ;
                    BMU(r+NREAD,:) = [x d_r2(x) start_r2(x) end_r2(x)] ;
%                 elseif (d1<d2 && d1<d3)
%                     if ( start1_r2>2360 && start1_r2<2501 && verbose )
%                         fprintf(1,'weight1 read1\n');
%                         printmat(weight{1},2500,2520) ;
%                         printmat(reads1(r).seqvect,2500-start1_r2,2500-start1_r2+20) ;
%                         printmat(seqvect1,2500,2520) ;
%                     end
%                     weight{1} = seqvect1 ;
%                     BMU(r,:) = [1 d1 start1_r1 end1_r1] ;
%                     BMU(r+NREAD,:) = [1 d1 start1_r2 end1_r2] ;
%                 elseif (d2<d3)
%                     if ( start2_r2>2360 && start2_r2<2501 && verbose )
%                         fprintf(1,'weight2 read1\n');
%                         printmat(weight{2},2500,2520) ;
%                         printmat(reads1(r).seqvect,2500-start2_r2,2500-start2_r2+20) ;
%                         printmat(seqvect2,2500,2520) ;
%                     end
%                     weight{2} = seqvect2 ;
%                     BMU(r,:) = [1 d2 start2_r1 end2_r1] ;
%                     BMU(r+NREAD,:) = [2 d2 start2_r2 end2_r2] ;
%                 else
%                     if ( start3_r2>2360 && start3_r2<2501 && verbose )
%                         fprintf(1,'weight3 read1\n');
%                         printmat(weight{3},2500,2520) ;
%                         printmat(reads1(r).seqvect,2500-start3_r2,2500-start3_r2+20) ;
%                         printmat(seqvect3,2500,2520) ;
%                     end
%                     weight{3} = seqvect3 ;
%                     BMU(r,:) = [1 d3 start3_r1 end3_r1] ;
%                     BMU(r+NREAD,:) = [3 d3 start3_r2 end3_r2] ;
                else
                    [~,x]=min(d_r2);
                    weight{x} = seqvect{x} ;
                    BMU(r,:) = [x d_r1(x) start_r1(x) end_r1(x)] ;
                    BMU(r+NREAD,:) = [x d_r2(x) start_r2(x) end_r2(x)] ;
                end
%             elseif (d1<d2 && d1<d3)
%                 if ( start1_r1>2360 && start1_r1<2501 && verbose )
%                     fprintf(1,'weight1 read1\n');
%                     printmat(weight{1},2500,2520) ;
%                     printmat(reads1(r).seqvect,2500-start1_r1,2500-start1_r1+20) ;
%                     printmat(seqvect1,2500,2520) ;
%                 end
%                 weight{1} = seqvect1 ;
%                 BMU(r,:) = [1 d1 start1_r1 end1_r1] ;
%                 [ d1, seqvect1, end1_r2, start1_r2 ] = DTWaverage( weight{1}, reads2(r).seqvect, 1, w2, 0, 1 ) ;
%                 BMU(r+NREAD,:) = [1 d1 start1_r2 end1_r2] ;
%                 if ( start1_r2>2360 && start1_r2<2501 && verbose )
%                     fprintf(1,'weight1 read2\n');
%                     printmat(weight{1},2500,2520) ;
%                     printmat(reads2(r).seqvect,2500-start1_r2,2500-start1_r2+20) ;
%                     printmat(seqvect1,2500,2520) ;
%                 end
%                 weight{1} = seqvect1 ;
%             elseif (d2<d3)
%                 if ( start2_r1>2360 && start2_r1<2501 && verbose )
%                     fprintf(1,'weight2 read1\n');
%                     printmat(weight{2},2500,2520) ;
%                     printmat(reads1(r).seqvect,2500-start2_r1,2500-start2_r1+20) ;
%                     printmat(seqvect2,2500,2520) ;
%                 end
%                 weight{2} = seqvect2 ;
%                 BMU(r,:) = [2 d2 start2_r1 end2_r1] ;
%                 [ d2, seqvect2, end2_r2, start2_r2 ] = DTWaverage( weight{2}, reads2(r).seqvect, 1, w2, 0, 1 ) ;
%                 BMU(r+NREAD,:) = [2 d2 start2_r2 end2_r2] ;
%                 if ( start2_r2>2360 && start2_r2<2501 && verbose )
%                     fprintf(1,'weight2 read2\n');
%                     printmat(weight{2},2500,2520) ;
%                     printmat(reads2(r).seqvect,2500-start2_r2,2500-start2_r2+20) ;
%                     printmat(seqvect2,2500,2520) ;
%                 end
%                 weight{2} = seqvect2 ;
%             else
%                 if ( start3_r1>2360 && start3_r1<2501 && verbose )
%                     fprintf(1,'weight3 read1\n');
%                     printmat(weight{3},2500,2520) ;
%                     printmat(reads1(r).seqvect,2500-start3_r1,2500-start3_r1+20) ;
%                     printmat(seqvect3,2500,2520) ;
%                 end
%                 weight{3} = seqvect3 ;
%                 BMU(r,:) = [3 d3 start3_r1 end3_r1] ;
%                 [ d3, seqvect3, end3_r2, start3_r2 ] = DTWaverage( weight{3}, reads2(r).seqvect, 1, w2, 0, 1 ) ;
%                 BMU(r+NREAD,:) = [3 d3 start3_r2 end3_r2] ;
%                 if ( start3_r2>2360 && start3_r2<2501 && verbose )
%                     fprintf(1,'weight3 read2\n');
%                     printmat(weight{3},2500,2520) ;
%                     printmat(reads2(r).seqvect,2500-start3_r2,2500-start3_r2+20) ;
%                     printmat(seqvect3,2500,2520) ;
%                 end
%                 weight{3} = seqvect3 ;
            else
                [~,x]=min(d_r1);
                weight{x} = seqvect{x} ;
                BMU(r,:) = [x d_r1(x) start_r1(x) end_r1(x)] ;
                if pe
                    [ d_r2(x), seqvect{x}, end_r2(x), start_r2(x) ] = DTWaverage( weight{x}, reads2(r).seqvect, 1, w2, 0, 1, cw, indel_cost, indel_weight ) ;
                    BMU(r+NREAD,:) = [x d_r2(x) start_r2(x) end_r2(x)] ;
                end
                weight{x} = seqvect{x} ;
            end
        end

        % Convert weight to nucleotide sequences
        weight_seq=struct('Header','','Sequence','','seqvect',0);
        for n=1:numw
            weight_seq(n).Header = ['weight' num2str(n)] ;
            weight_seq(n).seqvect = weight{n} ;
            weight_seq(n).Sequence = mat2nucleo(weight{n}) ;
            weight_seq(n).coverage = windowsa(weight{n}, 50, 0.99, BMU(BMU(:,1)==n,:), 0) ;
        end

        % Printing
        if 0
            ref.Sequence=reference.Sequence;ref.Header='weightZERO';ref.seqvect=reference.seqvect;
            ref0=ref_amplicons(1); ref0.Header='NC_027424'; ref0.seqvect=[.25; .25; .25; .25; 1];
            %seqalignviewer(multialign([weight2fasta(tree,weight) ref ref0]));
            tmp=rmfield(weight_seq,'coverage');
            true = fastaread('/mnt/bcd5a556-5a40-44d6-a73a-33acc92298e5/Pooling3_reps/geneious.fasta','TrimHeaders','true') ;
            seqalignviewer(multialign([tmp true(1:3) ref ref0]));
            fastawrite('./S10_reconstructed_rawreads_indels_pe.fa',weight_seq);
            %[~,ali] = nwalign(ref,ref0); showalignment(ali);
            %weight2fasta(tree,weight,'./S10_reconstructed_rawreads_indels.fa');
        end
    else % only one sequence to return ==> the root
        weight_seq=struct('Header','Nucleoveq_sequence','Sequence',reference.Sequence,'seqvect',reference.seqvect,'coverage',[]); % include "coverage" to keep consistent output
    end

    % % Window based analysis to get coverage
    % siz=200; ovlap=0.9; plotting=0;
    % for n=1:setsize
    %     [weight_seq(n).coverage,weight_seq(n).entropy] = windowsa(weight{n}, siz, ovlap, BMU(BMU(:,1)==n,:), plotting) ;
    % end

    if 0
        figure; plot(1:numel(wcovera1),[wcovera1; wcovera2; wcovera3],'-o'); legend({'weight1','weight2','weight3'});
        title(['window size ' num2str(siz) ', overlap ' num2str(100*ovlap) '%']);
        set(gca,'Xtick',1:(length(wcovera1)/10):length(wcovera1),'XTickLabel',xtl(round(1:(length(xtl)/10):end))) ;
        hold on; plot(1:numel(wcovera1),ones(numel(wcovera1)).*0.125*size(BMU,1)/(length(wcovera1)*(1-ovlap)),'red'); 
        text(numel(wcovera1),0.125*size(BMU,1)/(length(wcovera1)*(1-ovlap)),'12.5%'); hold off;
        hold on; plot(1:numel(wcovera1),ones(1,numel(wcovera1)).*0.250*size(BMU,1)/(length(wcovera1)*(1-ovlap)),'red'); 
        text(numel(wcovera1),0.250*size(BMU,1)/(length(wcovera1)*(1-ovlap)),'25%');hold off;
        hold on; plot(1:numel(wcovera1),ones(1,numel(wcovera1)).*0.625*size(BMU,1)/(length(wcovera1)*(1-ovlap)),'red'); 
        text(numel(wcovera1),0.625*size(BMU,1)/(length(wcovera1)*(1-ovlap)),'62.5%');hold off;
        hold on; plot(1:numel(wcovera1),ones(1,numel(wcovera1)).*1.000*size(BMU,1)/(length(wcovera1)*(1-ovlap)),'black'); hold off;
    end
    
    fprintf(1,'%s - Nucleoveq completed.\n', datestr(now)) ;

end
