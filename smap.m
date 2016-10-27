function [weight, BMU] = smap(nnode,epoch,LR,NS,reads,reference)
% SOM using DTWaverage

    % output
    filename = ['dtwaven_' num2str(epoch) '_' num2str(randi(1e6)) '.out'] ;
    fileid = fopen(filename,'w') ;
    
    % initialise map
    for n=1:nnode*nnode
        weight{n} = [ reference.seqvect(1:4,:); ones(1,size(reference.seqvect,2)) ] ;
        mutationpos = datasample(1:length(reference.seqvect(1:4,:)),length(reference.seqvect(1:4,:))*0.2,'Replace',false) ;
        for m=1:numel(mutationpos)
            weight{n}(1:4,m) = datasample(weight{n}(1:4,m),4,'replace',false) ;
        end
        entrop{n} = 0 ;
        wentrop{n} = 0 ;
        weightvarcov{n} = 0 ;
        nbmu(n) = 0 ; % record the cumulative number of hits
        nbmue(n) = 0 ; % record the number of hits for the current epoch
    end
    
    p = zeros(1,epoch) ; % mapping precision
    rposition = zeros(numel(reads),1) ;
    readlen = abs(mean(arrayfun(@(x) length(x.seqvect), reads))) ;
    BMU = zeros(numel(reads),2) ;
    
    % quantify map nodes
    fprintf(fileid,'#%g sequence matrices loaded\n',numel(reads)) ;
    %fprintf(fileid,'#Parameters -- %g, %g, %g, %g, %g, %g, %g, %g\n', epoch,LR(1),LR(2),NS(1),NC(1),NC(2),thexpan,gama ) ;
    for e=1:epoch
        fprintf(fileid,'\nepoch %s/%g %s -- ',sprintf(['%0' num2str(length(num2str(epoch))) 'd'],e),epoch,char(datetime('now')));
        % randomize the order which syllables are picked up
        numr=size(reads,2);
        rindex = randperm(numr) ;
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
        for s=rindex
            %%%%% restraint due to distance from BMU and learning restraint due to time
            % time constant that brings the learning rate and the neighborhood size to a very small value with iterations
            ED = exp(-e^2/(epoch*0.75)^2) ;
            %%%%% LEARNING RATE
            Ftime = max( LR(1)*ED , LR(2) ) ; % exponentially decreases from Ftime0
            % find BMU
            dist = +Inf ;
            for n=1:nnode*nnode
                [ d, ~, aligned_pos ] = DTWaverage( weight{n}, reads(s).seqvect, 1, 0.5, 0, 1 ) ;
                if (d<dist)
                    rposition(s) = aligned_pos ;
                    BMU(s,:) = [n d] ;
                    dist=d ;
                end
            end
            % mapping precision
            P = P + BMU(s,2) ;
            % update nodes
            bmu_coord = [ceil(BMU(s,1)/nnode) mod(BMU(s,1),nnode)] ;
            for L=1:nnode*nnode
                Edist = pdist([ ceil(L/nnode) mod(L,nnode) ; bmu_coord ]) ; % map structure
                Fdist = exp( (-Edist^2)/(2*(NeSt)^2) ) ;
                H = Ftime*Fdist ; % Learning Rate x Topological  Distance
                H = 1 - H ;
                %%%%% LEARNING STEP
                [ ~, updated_weight ] = DTWaverage( weight{L}, reads(s).seqvect, 1, H, 0, 1 ) ;
                weight{L} = updated_weight;
            end
            nbmu(BMU(s,1)) = nbmu(BMU(s,1))+1 ;
            nbmue(BMU(s,1)) = nbmue(BMU(s,1))+1 ;
        end
        p(e) = P/numel(reads) ;
        % test convergence
        converge = 1;
        se(e)=0;
        for w=1:nnode*nnode % all tips
            defined = round(0.95*size(weight{BMU(s,1)},2)) ; % what is shannon entropy if 95% of the sites would be fully defined
            undefined = size(weight{BMU(s,1)},2)-defined ;
            seqvectnull = [ repmat([1; 0; 0; 0],1,defined) repmat([.25; .25; .25; .25],1,undefined) ] ;
            senull = shannonEntropy(seqvectnull) ;
            sew = shannonEntropy(weight{w}) ;
            se(e) = se(e)+sew ;
            if sew>senull
                converge = 0;
            end
        end
        se(e) = se(e)/(nnode*nnode) ;
        if converge==1, fprintf(fileid,'\n#*** All weights are at minimum entropy95 ***\n'); break; end
    end
    
end