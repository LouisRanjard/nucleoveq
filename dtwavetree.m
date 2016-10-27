   classdef dtwavetree
   % write a description of the class here.

       properties
            root ;
            istrained ;
            epoch ;
            learningrate ;
            nstrength ;
            numchild ;
            thexpan ;
            gama ;
            fieldnam ;
            dico ;
            verbose ;
            ce ;
            fe ;
            reads ;
            mapprec ;
            size ;
            numleaves ;
            AIC = 2*numleaves-2*log(1/p(e)) ;
            BIC = -2*log(1/p(e))+numleaves*log(numsyll) ;
       end

       methods

           function obj = dtwavetree()
           % class constructor
               if(nargin > 0)
                    obj.root = node(0) ;
                    obj.istrained = false ;
                    obj.epoch = 1 ;
                    obj.learningrate = [0.99 0.99] ;
                    obj.nstrength = 0.4 ;
                    obj.numchild = [2 2] ;
                    obj.thexpan = 1e3 ;
                    obj.gama = 0.95 ;
                    obj.verbose = 0 ;
                    obj.ce = 0 ;
                    obj.fe = 1 ;
                    obj.reads = reads() ;
                    obj.mapprec ;
                    obj.size = 1 ;
                    obj.numleaves = 1 ;
                    obj.AIC = 0 ;
                    obj.BIC = 0 ;
               end
           end

           function train(obj)
           end

           function init_root(obj)
           end
           
           function load_reads(obj)
           end

           function get_mapprec(obj)
           end
           
           function get_AIC(obj)
           end

           function get_BIC(obj)
           end
           
           function plot(obj)
                figure;
                treeplot(tree) ;
                [x,y] = treelayout(tree);
                text(x,y,num2str([(1:length(tree))' ...
                nbmu' ...
                cellfun(@(x) shannonEntropy(x),weight)'...
                %arrayfun( @(x,y) siblings=sum(BMU(:,1)==find(tree==tree(y))); tipn=sum(BMU(:,1)==x); ,setdiff(1:numel(tree),unique(tree)),1:length(tree))...
                ],'%d\n%d\n%.3f'));
           end
           
       end
   end