function [tree, weight, BMU, features, nbmue, weightvarcov] = ETDTWrec(syllab,epoch,LR,NS,NC,thexpan,gama,fieldnam,dico,verbose,ce,fe,reference)
% performs a Evolving tree clustering analysis (Ranjard2008, Samarasinghe2006, Pakkanen2004)

% syllab.'fieldnam' are the matrices to classify
% epoch is the number of Epochs
% LR(1) is the initial Learning Rate decreasing until LR(2)
% NS is the neighborhood strength, if more than one element then use a
%     bell-shaped function instead of a constant value
%     NS(1) is the maximum value obtained at middle point of the learning process
%     NS(2) is the standard deviation of the function, if ==0 then NS(2) takes (#epochs)/5
% weight contains the weigth matrix of each index of tree
% NC(1) is the initial number of children per node and NC(2) is the final number of children
% thexpan is the threshold above which a node is splitted
%       only nodes which have a cluster scatter above the mean can be divided
% gama is the weight decay on the hit counter (Pakkanen2006)
% fielnam is the name of the field in syllab
% dico is a classification to be used to improve the convergence of the clustering
%       must be decreasingly ordered according to the distance
% verbose: verbose mode? 1 or 0
% example ETDTWrec(syllab,5,[0.9 0.01],[4 1],[2 2],numel(syllab)*0.2,0.95,'cepvect',[],1)

% tree contains for each index the index of its father
% weight contains the weight matrices (cluster centres)
% BMU records the path for each sample in the ETree

% LEARNING AND DISPLAY PARAMETERS
if nargin<13
    reference = [] ;
    if nargin<12
        fe = 0 ;
        if nargin<11
            ce = 0 ;
            if nargin<10
                verbose = 0 ;
                if nargin<9
                    dico = [] ;
                    if nargin<8
                        fieldnam = 'cepvect' ;
                        if nargin<7
                            gama = 1 ;
                            if nargin<6
                                thexpan = round(numel(syllab)/10) ; % arbitrarily set the splitting threshold
                                if nargin<5
                                    NC = [2 2] ; % default is binary trees
                                    if nargin<4
                                        NS = [3 3] ;
                                        if nargin<3
                                            LR = [0.9 0.01] ; % minimal learning rate (Samarasinghe2006 uses a threshold of 0.01)
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if numel(reference)>0 && ~isfield(reference,'varcov')
    error('Field "varcov" not found in reference');
end

% debugging mode
DEBUG=0;

% output
rng('shuffle'); % creates a different seed each time
filename = ['dtwaven_' num2str(epoch) '_' num2str(randi(1e6)) '.out'] ;
fileid = fopen(filename,'w') ;
OUTFILE=0;
if OUTFILE
    filelog=fopen([filename '.log'],'w') ;
    filelog2=fopen([filename '_2.log'],'w') ;
end

% to make computation faster, record the dico index of each syllable
if numel(dico)>0
    syldico=zeros(1,numel(syllab)) ;
    for si=1:numel(syllab)
        for di=1:numel(dico)
            if find(dico(di).example==si)
                syldico(si)=di ;
                break ;
            end
        end
    end
end


% INITIALISATION
%%%%% mean sizes of syllab cepstrum to find the size of weight matrices
% max 100 random samples
idxs = ceil(numel(syllab).*rand(1,min(100,numel(syllab)))) ; % avoid extrem values
% fprintf(1,'%d %d\n',size([syllab(idxs).(fieldnam)],1),size([syllab(idxs).(fieldnam)],2)) ;
vals = reshape([syllab(idxs).(fieldnam)],1,[]) ;
x = size([syllab.(fieldnam)],1) ; % frequency mean size
if (fe)
    y = fe ;
else     
    y = round(size([syllab.(fieldnam)],2)/numel(syllab)) ; % time mean size
end

% Tree initialisation with random weights values picked from vals
tree(1) = 0 ; % Root node
numchildren = NC(1) ;
if fe && numel(reference)>0
    %weight{1} = reference ;
    weight{1} = reference.seqvect(1:4,:) ;
    weightvarcov{1} = reference.varcov ;
    entrop{1} = shannonEntropy(weight{1}) ;
    wentrop{1} = shannonEntropy(weight{1})/numel(syllab) ;
end
nbmu(1) = numel(syllab) ; % record the cumulative number of hits
nbmue(1) = numel(syllab) ; % record the number of hits for the current epoch
for w=1:numchildren
    tree(w+1) = 1 ; % its father is the root
    node(w+1).parent = 1 ;
    rindex = 1 + round((size(vals,2)-1)*rand(x,y)) ;
    nbmu(w+1) = 0 ; % will store the number of time this node is BMU (hit counter)
    nbmue(w+1) = 0 ; % will store the number of time this node is BMU (hit counter) in the current epoch
    if (fe)
        if numel(reference)==0
            weight{w+1} = randmatseq(fe) ; % for nucleotide sequences, random sequence
        else
            %weight{w+1} = reference ; % initialise the tree with the input reference sequence
            weight{w+1} = [ weight{1}; ones(1,size(reference.seqvect,2)) ] ; % add one line for persistence of a position
        end
    else
        weight{w+1} = vals(rindex) ; % +1 because weight(1) is the root
    end    
    entrop{w+1} = shannonEntropy(weight{w+1}) ;
    wentrop{w+1} = shannonEntropy(weight{w+1})/(numel(syllab)/2) ;
    weightvarcov{w+1} = 0 ;
end
clear vals ;

% DISPLAY?
if verbose==1
    %%% figures
    f = figure('Units','pixels','Position',[300 100 600 840]) ;
    Panel1 = uipanel('Units','pixels','Position',[20 20 270 180],'Title','Davies-Bouldin','Parent',f) ;
    Panel2 = uipanel('Units','pixels','Position',[310 20 270 180],'Title','Evolving Tree','Parent',f) ;
    Panel3 = uipanel('Units','pixels','Position',[20 220 270 180],'Title','Furthest Sample Mean Distance','Parent',f) ;
    Panel4 = uipanel('Units','pixels','Position',[310 220 270 180],'Title','Mapping Precision','Parent',f) ;
    Panel5 = uipanel('Units','pixels','Position',[20 420 270 180],'Title','Closest Sample Mean Distance','Parent',f) ;
    Panel6 = uipanel('Units','pixels','Position',[310 420 270 180],'Title','Inner Scatter Distribution','Parent',f) ;
    Panel7 = uipanel('Units','pixels','Position',[20 620 270 180],'Title','Dividing Events','Parent',f) ;
    mapprec = axes('parent',Panel4) ;
    treep = axes('parent',Panel2) ;
    daviesbouldinax = axes('parent',Panel1) ;
    innerscatterax = axes('parent',Panel6) ;
    closestsample = axes('parent',Panel5) ;
    dividingevents = axes('parent',Panel7) ;
    maxclustdist = axes('parent',Panel3) ;
elseif usejava('desktop') && false
    % shows progress
    h = waitbar(0,'learning...') ;
end

% print the parameter values
fprintf(fileid,'#%g sequence matrices loaded\n',numel(syllab)) ;
if isfield(reference,'varcov')
    fprintf(fileid,'#reference.varcov %.4f\n',reference.varcov) ;
end
if isfield(reference,'entropy')
    fprintf(fileid,'#reference.entropy %.4f\n',reference.entropy) ;
end
if numel(NS)>1
    fprintf(fileid,'#Learning parameters %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n', epoch,LR(1),LR(2),NS(1),NS(2),NC(1),NC(2),thexpan,gama ) ;
else
    fprintf(fileid,'#Learning parameters %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n', epoch,LR(1),LR(2),NS(1),NC(1),NC(2),thexpan,gama ) ;
end
if verbose==1
    fprintf(fileid,'#Epoch\t Time\t Neighborhood Size\t Mapping precision\t Davies-Bouldin\t Tree size\t Leaves number\n' ) ;
else
    fprintf(fileid,'Epoch Time Neighborhood_Size Mapping_precision Mean_leaf_entropy Joint_leaf_entropy Tree_size Leaves_number' ) ;
end

% some memory allocation
if DEBUG
    BMU = zeros(numel(syllab),1,2) ; % Path to Best Matching Unit and Score, can't know the max depth yet so just initialize with a value of 1
else
    BMU = zeros(numel(syllab),2) ; % Best Matching Unit and Score only
end
p = zeros(1,epoch) ; % mapping precision
DB = [] ; % Davies Bouldin index
meanCS = [] ; % record the mean distance of the closest sample from each cluster centre (leaf)
meanFS = [] ; % record the mean distance of the furthest sample from each cluster centre
dividing = [] ; % record when there is a dividing event during learning
divided = 0 ; % test to know when there has been dividing event
%features = cell2table(cell(epoch,5)) ; % store features of the Etree
%features.Properties.VariableNames = {'TreeSize' 'NumLeaves' 'MappingPrecision' 'AIC' 'BIC'} ;
features = struct('TreeSize',zeros(1,epoch),'NumLeaves',zeros(1,epoch),'MappingPrecision',zeros(1,epoch),...
    'AIC',zeros(1,epoch),'BIC',zeros(1,epoch),'AllEntropy',zeros(1,epoch)) ;
rposition = zeros(numel(syllab),1) ;
readlen = abs(mean(arrayfun(@(x) length(x.seqvect), syllab))) ;

% get total number of steps to implement a step counter
steprecord = ceil(epoch*numel(syllab)*0.01) ;
stepcount = 0 ;

% use a bell shaped neighborhood size?
if numel(NS)>1
    if NS(2)==0 % no standard deviation input, take the (total number of epoch)/5
        NS(2) = epoch/5 ;
    end
    % compute the factor to be used during updates
    nsfactor = NS(1)/(1/(NS(2)*sqrt(2*pi))) ;
end    

for e=1:epoch % each epoch
    if 0 % some plotting
        figure;
        treeplot(tree) ;
        [x,y] = treelayout(tree);
        text(x,y,num2str([(1:length(tree))' ...
                          nbmu' ...
                          cellfun(@(x) shannonEntropy(x),weight)'...
                          %arrayfun( @(x,y) siblings=sum(BMU(:,1)==find(tree==tree(y))); tipn=sum(BMU(:,1)==x); ,setdiff(1:numel(tree),unique(tree)),1:length(tree))...
                          ],'%d\n%d\n%.3f'));
        drawnow;
        %pause;
    end
    %currtime=fix(clock);
    %fprintf(1,'\nepoch %g/%g %g-%g-%g %g:%g -- ',e,epoch,currtime(3),currtime(2),currtime(1),currtime(4),currtime(5));
    % randomize the order which syllables are picked up
    numsyll=size(syllab,2);
    sindex = randperm(numsyll) ;
    %[~,sindex] = sort(BMU(:,dep,2),'descend') ; % to be tested
    P = 0 ;
    %%%% compute the neighboring size for this epoch
    if numel(NS)>1
        NeSt = max( 1 , (1/(NS(2)*sqrt(2*pi)))*exp(-(((e)-((epoch+1)/2)).^2)/(2*(NS(2))^2))*nsfactor ) ;
        NeSt = round(NeSt) ;
    else
        NeSt = NS ;
    end
    % reset the number of hits per epoch
    nbmue = nbmue .* 0 ;
    % every syllable
    for s=sindex
        if DEBUG && fe
            for L=2:numel(tree)
                if sum(sum(weight{L},1)~=1)>0
                    1;
                end
            end
        end
        %%%%% restraint due to distance from BMU and learning restraint due to time
        % time constant that brings the learning rate and the neighborhood size to a very small value with iterations
        ED = exp(-e^2/(epoch*0.75)^2) ;
        %%%%% LEARNING RATE
        Ftime = max( LR(1)*ED , LR(2) ) ; % exponentially decreases from Ftime0
        %%%%% if free-end (alignment to reference), parse true position
        if DEBUG && fe
            original_pos = str2double(syllab(s).Header(strfind(syllab(s).Header,'pos')+3:end)) ;
        end
        %%%%% search for BMU in the tree
        if DEBUG % store the full path to BMU
            tmp = findBMU(syllab(s).seqvect, tree, weight, ce, fe, 'full') ;
            BMU(s,:,:) = 0 ;  % re-initialised BMU(s,:,:) if the new BMU in the tree is not as deep as the previous one
            BMU(s,1:size(tmp,1),1:2) = tmp ;
            %%%%% mapping precision
            P = P + BMU(s,dep,2) ;
        else % only store last BMU index and score
            [ treepath, aligned_pos ] =  findBMU(syllab(s).seqvect, tree, weight, ce, fe);
            BMU(s,:) = treepath ;
            rposition(s) = aligned_pos ;
            %%%%% mapping precision
            P = P + BMU(s,2) ;
        end
        %%%%% Use the samples in the same cluster to update (according to the dico)
%         if numel(dico)>0 && syldico(s)>0
%             cws = dico(syldico(s)).example ;
%             for cwi=1:numel(dico(syldico(s)).example)
%                 for L=2:numel(tree)
%                     %%%%% NEIGHBORHOOD SIZE
%                     Edist = TreeDist(tree,BMU(s,dep,1),L) ;
%                     %Fdist = exp( (-Edist^2)/(2*(NeSt*exp(-e/(epoch*0.5)))^2) ) ; % old bug
%                     Fdist = exp( (-Edist^2)/(2*(NeSt)^2) ) ;
%                     H = Ftime*Fdist ;
%                     %%%%% LEARNING STEP
%                     %if ( H>0.01 ) % in order to save time if the factor is too small do not compute
%                         [ ~, weight{L} ] = DTWaverage( weight{L}, syllab(cws(cwi)).(fieldnam), 1, H, ce, fe ) ;
%                     %end
%                 end
%             end
%         end
        %%%%% UPDATE the tree except the root in the neighbourhood of the Best Matching Unit by pulling them closer to the input matrix
        % H is the weight for the first matrix sent to DTWaverage (weight{L}), 
        % hence it needs to increase with distance from BMU (Edist) ???BUG???
        for L=2:numel(tree)
                %%%%% NEIGHBORHOOD SIZE
                if DEBUG
                    Edist = TreeDist(tree,BMU(s,dep,1),L) ;
                else
                    Edist = TreeDist(tree,BMU(s,1),L) ;
                end
                %Fdist = exp( (-Edist^2)/(2*(NeSt*exp(-e/(epoch*0.5)))^2) ) ; % old bug
                Fdist = exp( (-Edist^2)/(2*(NeSt)^2) ) ;
                H = Ftime*Fdist ; % Learning Rate x Tree  Distance
                if fe
                    H = 1 - H ;
                end
                %%%%% LEARNING STEP
                %if ( H>0.01 ) % in order to save time if the factor is too small do not compute
                    [ ~, updated_weight, aligned_pos ] = DTWaverage( weight{L}, syllab(s).(fieldnam), 1, H, ce, fe ) ;
                    if DEBUG && fe && L==BMU(s,dep,1) && ( aligned_pos~=(original_pos+size(syllab(s).(fieldnam),2)-1) ) % aligned_pos is index in C (start at 0)
                        %mat2nucleo( syllab(s).(fieldnam) )
                        %mat2nucleo( weight{L}(:,original_pos:(original_pos+size(syllab(s).(fieldnam),2)-1)) )
                        %mat2nucleo( weight{L}(:,(aligned_pos-size(syllab(s).(fieldnam),2)+1):aligned_pos) )
                        %mat2nucleo( updated_weight(:,original_pos:(original_pos+size(syllab(s).(fieldnam),2)-1)) )
                    end
                    weight{L}=updated_weight;
                %end
        end
        %%%%% increases nbmu 
        if DEBUG
            nbmu(BMU(s,dep,1)) = nbmu(BMU(s,dep,1))+1 ;
        else
            nbmu(BMU(s,1)) = nbmu(BMU(s,1))+1 ;
            nbmue(BMU(s,1)) = nbmue(BMU(s,1))+1 ;
        end
        %%%%% Update weighted entropy of BMU (keep a window of the last values only)
        %swindow = 100 ;
%         swindow = nbmu(tree(BMU(s,1)))/3 ;
%         if length(entrop{BMU(s,1)})==swindow
%             entrop{BMU(s,1)}(1:(swindow-1)) = entrop{BMU(s,1)}(2:swindow) ;
%             entrop{BMU(s,1)}(swindow) = shannonEntropy(weight{BMU(s,1)}) ;
%             wentrop{BMU(s,1)}(1:(swindow-1)) = entrop{BMU(s,1)}(2:swindow) ;
%             wentrop{BMU(s,1)}(swindow) = shannonEntropy(weight{BMU(s,1)})/nbmu(BMU(s,1)) ;
%         else
%             entrop{BMU(s,1)}(end+1) = shannonEntropy(weight{BMU(s,1)}) ;
%             wentrop{BMU(s,1)}(end+1) = shannonEntropy(weight{BMU(s,1)})/nbmu(BMU(s,1)) ;
%         end
         
        %%%%% DIVIDING STEP
        % compute cutoff entropy
        % -> calculate the Shannon entropy of a sequence where 95% of the sites are fully defined
        wlength = size(weight{BMU(s,1)},2) ;
        defined = round(0.95*wlength) ; 
        undefined = size(weight{BMU(s,1)},2)-defined ;
        seqvectnull = [ repmat([1; 0; 0; 0],1,defined) repmat([.25; .25; .25; .25],1,undefined) ] ;
        cutoff_e = shannonEntropy(seqvectnull) ;
        % get statistics on coverage
        [ weightcov, varcov, ~ ] = wcoverage( weight{BMU(s,1)}, readlen, rposition(BMU(:,1)==BMU(s,1)) ) ;
        weightvarcov{BMU(s,1)} = [weightvarcov{BMU(s,1)} varcov] ;
        % set factors for criteria to subdivise a node
        nbmu_factor=2; % default 2
        entropy_factor=1; % default 1
        weightcov_factor=0.5; % default 0.50
        varcov_factor=10; % default 2
        %if nbmu(BMU(s,1))>thexpan                                                                                          % Use a hit counter
        %if shannonEntropy(weight{BMU(s,1)}) > shannonEntropy(weight{tree(BMU(s,1))})                                        % Entropy is greater than parent's
            %disp(cellfun(@(x) shannonEtree(x),weight)) ; disp(BMU(s,1)) ;
        %if ( length(entrop{BMU(s,1)})==swindow && entrop{BMU(s,1)}(end)>entrop{BMU(s,1)}(1) )                              % wEntropy is greater than previous entropy
        %if ( length(entrop{BMU(s,1)})==swindow && entrop{BMU(s,1)}(end)==max(entrop{BMU(s,1)}) )                           % wEntropy is greater than n previous values
        %if ( length(entrop{BMU(s,1)})==swindow && entrop{BMU(s,1)}(end)>(mean(entrop{BMU(s,1)})+2*var(entrop{BMU(s,1)})) ) % wEntropy is greater than ( n previous values + 2*variance )
        %if ( length(entrop{BMU(s,1)})==swindow &&...
        %        entrop{BMU(s,1)}(end)>(mean(entrop{BMU(s,1)})+2*var(entrop{BMU(s,1)})) &&...
        %        nbmu(BMU(s,1))>nbmu(tree(BMU(s,1))) )                                                                      % same + number of BMU of node is greater than parent 
        %if ( nbmu(BMU(s,1))>nbmu(tree(BMU(s,1))) )                                                                         % hit counter greater than parent's (NOT BAD)
        %if ( nbmu(BMU(s,1)) >= ( nbmu(tree(BMU(s,1))) / 2 )  && ...                                                          % hit counter greater than parent's /2
        %        entrop{BMU(s,1)}(end) > entrop{tree(BMU(s,1))}(end) )
        % AND wEntropy greater than parent's (never true)
        if OUTFILE
            fprintf(filelog,'%d\t%d\t%d\t%d\n',nbmu(BMU(s,1)),shannonEntropy(weight{BMU(s,1)}),sum(weightcov==0),varcov) ;
            fprintf(filelog2,'%d\t%d\t%d\t%d\n',nbmu(BMU(s,1))>=(nbmu(tree(BMU(s,1)))/nbmu_factor), shannonEntropy(weight{BMU(s,1)})>(entropy_factor*cutoff_e),...
                                        sum(weightcov==0)<(weightcov_factor*wlength), varcov<(varcov_factor*reference.varcov)) ;
        end
        if ( nbmu(BMU(s,1))>=(nbmu(tree(BMU(s,1)))/nbmu_factor) &&...
                shannonEntropy(weight{BMU(s,1)})>(entropy_factor*cutoff_e) &&...
                sum(weightcov==0)<(weightcov_factor*wlength) &&...
                varcov<(varcov_factor*reference.varcov) )
        %        entrop{BMU(s,1)}(end) > (mean(entrop{BMU(s,1)})+1*var(entrop{BMU(s,1)})) )                                 % AND wEntropy is greater than 2 z-score
        %        entrop{BMU(s,1)}(end) > ( entrop{tree(BMU(s,1))}(end) / 2 ) )                                              % AND wEntropy greater than parent's /2 (BAD : always true)
        % [1:numel(weight) ; nbmu ; cellfun(@(x) x(end),entrop) ; cellfun(@(x) x(end),wentrop)]
            for w=1:numchildren
                tree = [tree BMU(s,1)] ;
                if (fe)
                    weight{numel(tree)} = weight{BMU(s,1)} ; % simply copy parent node
%                     entrop{numel(tree)} = shannonEntropy(weight{numel(tree)}) ;
%                     wentrop{numel(tree)} = shannonEntropy(weight{numel(tree)})/(nbmu(tree(BMU(s,1)))/2) ;
                    weightvarcov{numel(tree)} = 0 ;
                end
                nbmu = [nbmu 0] ;
                nbmue = [nbmue 0] ;
            end
        %elseif nbmu(BMU(s,1)) >= ( nbmu(tree(BMU(s,1))) / 2 )
        %    disp(false);
        end
        %%%%% computes WITHIN CLUSTER INDEX (DaviesBouldin clustering index)
        %clustidx = unique(BMU(BMU(:,:,1)>0)) ; % the nodes used in the tree !!DO NOT WORK IN OCTAVE!!
%         BMU1 = BMU(:,:,1) ;
%         clustidx = unique(BMU1(BMU1>0)) ;
%         leaves = setdiff(1:numel(tree),unique(tree)) ; % the leaves i.e potentially clusters in the tree
%         clustidx = intersect(clustidx,leaves) ; % all the nodes which are also leaves
%         if numel(clustidx)>0
%             S = zeros(1,max(clustidx)) ; % within cluster scatter
% %             if verbose==1
% %                 Smaxtmp = zeros(1,max(clustidx)) ; % maximum cluster distance %noneed%
% %                 D = zeros(max(clustidx),max(clustidx)) ; % between cluster separation (for DAVIES-BOULDIN) %noneed%
% %             end
%             for n=1:numel(clustidx)
%                cidx = clustidx(n) ;
%                % within cluster scatter
%                [a, b] = find(BMU(:,:,1)==cidx) ;
%                S(cidx) = sum( BMU(a,b(1),2) ) / numel(a) ; % only consider b(1) because they're all the same depth
%                %%%%% DIVIDING STEP
%                % divide if S>mean(mapping precision) and expansion threshold reached (except for first epoch)
%                if e==1 && nbmu(cidx)>thexpan
%                    divid=1;
%                elseif  nbmu(cidx)>thexpan && S(cidx)>(P(P>0)/numel(P(P>0)))
%                    divid=1;
%                else
%                    divid=0;
%                end
%                if divid==1
%                     divided = 1 ; % record that there has been dividing for at least one cluster (for display)
%                     for w=1:numchildren
%                         tree = [tree cidx] ;
%                         if (fe)
%                             weight{numel(tree)} = weight{cidx} ; % simply copy parent node
%                         else
%                             weight{numel(tree)} = (rand(size(weight{cidx},1),size(weight{cidx},2))*0.20+0.9).*weight{cidx} ; % add some random noise
%                         end
%                         nbmu = [nbmu 0] ;
%                     end
%                     % initialise the new child
%                     [a, b] = find(BMU(:,:,1)==cidx) ;
%                     for m=1:numel(a)
%                        BMU(a(m),b(m)+1,1) = numel(tree)-floor(rand*numchildren) ; % randomly assign a child
%                        BMU(a(m),b(m)+1,2) = BMU(a(m),b(m),2) ; % use the same distance 
%                        % this makes DAVIES-BOULDIN worse because from one cluster with a bad inner scatter, 
%                        % two are created. this is specially true at the begining of the learning process when the
%                        % number of cluster is quite small.
%                     end
%                end
%                % between cluster separation (for DAVIES-BOULDIN) %noneed%
% %                if verbose==1
% %                    clustidxdif = clustidx(clustidx~=cidx) ;
% %                    for n2=1:numel(clustidxdif)
% %                        cidx2 = clustidxdif(n2) ;
% %                        D(cidx,cidx2) = DTWaverage( weight{cidx}, weight{cidx2} ) ;
% %                    end
% %                end
%             end
%             % DISPLAY?
%             if verbose==1
%                 if stepcount==steprecord
%                     stepcount = 0 ;
%                     % COMPUTE SOME CLUSTERING STATISTICS every 0.01% update steps
%                     % need to recompute all the distances because there has been updates (so not possible to use BMU)
%                     [distSN, distNN] = distSNN(syllab,weight,'',fieldnam,'SNN',0) ;
%                     [x,y,maxdep,s] = treelayout(tree) ;
%                     % depth of each node
%                     depth = zeros(1,numel(tree)) ;
%                     for n = 1:numel(tree)
%                         a=n ;
%                         while ( tree(a)~=0 ) depth(n)=depth(n)+1; a=tree(a); end
%                     end
%                     leaves = setdiff(1:numel(tree),unique(tree)) ;
%                     % the cluster centres are the ones of depth==limdep and the ones of inferior depth but which are leaves
%                     diconodes = [ find(depth==dep) intersect(find(depth<dep),leaves) ] ;
%                     % defines the cluster (dictionary)
%                     dico = struct('example',cell(1,numel(diconodes))) ;
%                     % redefine distSN and distNN with just the cluster centroids
%                     distSNc = distSN(:,[1 diconodes]) ;
%                     % build dico
%                     for s=1:size(distSNc,1)
%                         dicoidx = find(distSNc(s,2:end)==min(distSNc(s,2:end)),1) ; % add 1 because we avoided the root
%                         dico = AddExampleDico( dico, dicoidx, s ) ;
%                     end
%                     [ cs1, fs1 ] = ClustStat(dico,diconodes,distSN) ;
%                     db1 = DaviesBouldin2(dico,diconodes,distSN,distNN) ;
%                     %%% PLOT THESE STATISTICS
%                     meanCS = [meanCS cs1]; axes(closestsample); plot(meanCS);
%                     meanFS = [meanFS fs1]; axes(maxclustdist); plot(meanFS);
%                     DB = [DB db1]; axes(daviesbouldinax); plot(DB);
%                     dico = [] ;
%                 else
%                     stepcount = stepcount+1 ;
%                 end
% 
%                 % dividing events
%                 if divided==1 % there has been dividing
%                     divided=0;
%                     dividing = [dividing 1] ;
%                 else
%                     dividing = [dividing 0] ;
%                 end
%                 axes(dividingevents);
%                 plot(dividing); 
%                 drawnow;
% %                 % closest sample mean distance from each cluster centre %noneed%
% %                 distSN = distSNN(syllab,weight([1 clustidx]),'',fieldnam,'SN',0) ; % watch out, distSNN() skips the first weight which is supposly the root of the ETree
% %                 meanCS = [meanCS sum(min(distSN,[],1))/numel(clustidx)] ;
% %                 axes(closestsample);
% %                 plot(meanCS);
% %                 % furthest sample mean distance from each cluster centre (need to know the clusters) %noneed%
% %                 FS = zeros(1,numel(clustidx)) ;
% %                 for cidx3=1:numel(clustidx)
% %                     [a b] = find(BMU(:,:,1)==clustidx(cidx3)) ;
% %                     FS = max(distSN(a,1+cidx3)) ;
% %                 end
% %                 meanFS = [meanFS sum(FS)/numel(clustidx)] ;
% %                 axes(maxclustdist);
% %                 plot(meanFS); 
% %                 drawnow;
%                 % distribution of inner scatter
%                 axes(innerscatterax);
%                 hist(S(S>0),10);
% %                 % between cluster separation (for DAVIES-BOULDIN) %noneed%
% %                 R = zeros(1,max(clustidx)) ; % store the max ratio of S and D
% %                 for n=1:numel(clustidx)
% %                     cidx = clustidx(n) ;
% %                     indexdif = clustidx(clustidx~=cidx) ;
% %                     R(cidx) = max( (S(cidx)+S(indexdif))./D(cidx,indexdif) ) ;
% %                 end
% %                 % compute DAVIES-BOULDIN index %noneed%
% %                 DB(find(DB==0,1)) = sum(R)/numel(clustidx) ;
% %                 axes(daviesbouldinax); 
% %                 plot(DB(DB>0)); 
% %                 drawnow;
%             elseif usejava('desktop') && false
%                 waitbar((s+(e-1)*numel(syllab))/(numel(syllab)*epoch)) ;
%             end
%         end
%         %%%%%%%%% plot the tree with labels
        % DISPLAY?
        if verbose==1
            axes(treep) ; 
            treeplot(tree) ;
            [x,y] = treelayout(tree);
            %text(x,y,num2str((1:length(tree))')); % index
            %text(x,y,num2str(nbmu(ii)')); % number of bmu
            text(x,y,num2str([(1:length(tree))' nbmu'])); % both
            drawnow ;
        end
    end
    % weight decay on hit counter (Pakkanen2006) after each epoch
    nbmu = round(gama.*nbmu) ;
    % expansion size decay (Louis2007)
    numchildren = max( NC(1)-(round(NC(1)/epoch)) , NC(2) ) ;
    % Mapping precision
    tips = setdiff(1:numel(tree),unique(tree)) ;
    numtips = numel(tips) ;
    p(e) = P/numsyll ;
    % Shannon Entropy, control and stop if all nodes have low entropy
    converge = 0; %%% DEBUG %%% prevent convergence here, should be "converge = 1;"
    se(e)=0;
    sew = zeros(1,numel(weight)) ;
    for w=tips % all tips
        defined = round(0.95*size(weight{BMU(s,1)},2)) ; % what is shannon entropy if 95% of the sites would be fully defined
        undefined = size(weight{BMU(s,1)},2)-defined ;
        seqvectnull = [ repmat([1; 0; 0; 0],1,defined) repmat([.25; .25; .25; .25],1,undefined) ] ;
        senull = shannonEntropy(seqvectnull) ;
        sew(w) = shannonEntropy(weight{w}) ;
        se(e) = se(e)+sew(w) ;
        if sew(w)>senull
            converge = 0;
        end
    end
    se(e) = se(e)/numtips ;
    % Clustering criterion
    %fprintf(1,'BIC=%g, ', -2*log(1/p(e))+numtips*log(numsyll) );
    %fprintf(1,'AIC=%g, ', 2*numtips-2*log(1/p(e)) );
    %% Prune tree for zero hits tips
    %[tree,weight,nbmue] = pruneTip(tree,weight,nbmue,find(nbmue>0)) ; %%% DEBUG %%%
    %% Fortify the tree by realigning the reads to their BMU
    % choose pulling weight so that, if for a given position all the N reads r are different
    % from reference R, then the position in reference is converted to the read
    % position value after going through all the reads exactly. It cannot reach one so we use 0.9999
    % pw = (1-0.9999)^(1/N)
    pw = 0.0001^(1/numel(syllab)) ;
    %[weight, BMU] = fortify(10, reference.entropy, syllab, weight, pw, BMU, 0, 1, 1) ; % 10 iterations max unless entropy is reached
    [weight, BMU] = fortify(10, 0, syllab, weight, pw, BMU, 0, 1, 1) ; % 10 iterations max unless entropy is reached
    
    %% Save some features of the Etree
    features.TreeSize(e) = numel(tree) ;
    features.NumTips(e) = numtips ;
    features.MappingPrecision(e) = p(e) ;
    features.MeanShannonEntropy(e) = se(e) ;
%     features.AIC(e) = 2*numtips-2*log(1/p(e)) ;
%     features.BIC(e) = -2*log(1/p(e))+numtips*log(numsyll) ;
    features.AllEntropy(e) = shannonEntropy_s(weight(tips)) ;
    fprintf(fileid,'\n%s/%g\t %s\t ',sprintf(['%0' num2str(length(num2str(epoch))) 'd'],e),epoch,char(datetime('now'))); % format the printing of "e" to have the same width as "epoch"
    if verbose==1
        fprintf(fileid,'%f\t %f\t %f\t %g\t %g', NeSt, p(e), DB(end), numel(tree), numtips ) ;
    else
        fprintf(fileid,'%.2f\t %.4f\t %.4f\t %.4f\t %g\t %g', NeSt, p(e), se(e), features.AllEntropy(e), numel(tree), numtips ) ;        
    end
    
    %% check some criteria to stop learning
    % Stop if all nodes have low entropy (mean entropy)
    if converge==1, fprintf(fileid,'\n#*** All weights are at minimum entropy95 ***\n'); break; end
    % Check overall entropy versus the initial entropy of the root
    if (e>2 && features.AllEntropy(e)>reference.entropy), fprintf(fileid,'\n#*** Joint tips entropy is greater than root entropy ***\n'); break; end
    % Check if overall entropy is decreasing (=overfitting)
    %if (e>2 && features.AllEntropy(e)<features.AllEntropy(e-2)),  fprintf(fileid,'\n#*** Joint tips entropy is decreasing ***\n'); break; end
    % Stopping criteria NOT met -> Force the dividing of the highest entropy node
    [~,hottip] = max(sew) ;
    for w=1:numchildren
        tree = [tree hottip] ;
        if (fe)
            weight{numel(tree)} = weight{hottip} ; % simply copy parent node
            weightvarcov{numel(tree)} = 0 ;
        end
        nbmu = [nbmu 0] ;
        nbmue = [nbmue 0] ;
    end
    % DISPLAY?
    if verbose==1, axes(mapprec); plot(p); drawnow; end
    %%%%%%%%%%%%%%%% SAVING THE TEMPORARY TREE %%%%%%%%%%%%%%%%
    %save( '-v7', './tree_tmp.mat', 'tree', 'weight','e','nbmu' );
end

%if verbose==0 && usejava('desktop')
%    close(h) ;
%end
%delete('./tree_tmp.mat') ;
if verbose==1 && usejava('desktop') % plot the tree with nBMU
    figure; 
    treeplot(tree) ;
    [x,y] = treelayout(tree);
    %text(x,y,num2str([(1:length(tree))' nbmu']));
    text(x,y,num2str([nbmu']));
end
fprintf(fileid,'\n#------------------------------------------------\n');
fclose(fileid);
if OUTFILE
    fclose(filelog);
    fclose(filelog2);
end

end
