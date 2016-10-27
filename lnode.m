classdef lnode < handle
% write a description of the class here.

   properties
       idl ;
       parent ;
       child ;
       weight ;
       entropy ;
       nbmu ;
       mappingprecision ;
       coverage ;
   end

   methods

       % CONSTRUCTOR
       function obj = lnode(parent,weight,identifier)
       % class constructor
           if(nargin==0)
             obj.parent = 0 ;
             obj.child = [0 0] ;
             obj.weight = {} ;
             obj.entropy = 0 ;
             obj.nbmu = 0 ;
             obj.mappingprecision = 0 ;
           else
             obj.parent = parent ;
             obj.weight = weight ;
             obj.idl = identifier ;
             obj.weight.Header = ['sequence_' num2str(obj.idl)] ;
             obj.weight.varcov = nan ;
             obj.child = nan ;
             obj.entropy = nan ;
             obj.nbmu = nan ;
             obj.mappingprecision = nan ;
           end
       end

       %function align(obj,read)
           % DTWaverage()
       %end
       
       function obj = update(obj,read,factor)
           %se = shannonEntropy(obj.weight(:,read.position:(read.position+size(read.seqvect,2)))) ;
           for position = read.position:(read.position+size(read.seqvect,2)-1)
               if max(obj.weight.seqvect(1:4,position))<0.9
                   obj.weight.seqvect(1:4,position) = read.seqvect(:,position-read.position+1).*factor + obj.weight.seqvect(1:4,position).*(1-factor) ;
               end
           end
       end
       
       function obj = accord(obj,reads)
           obj.entropy = shannonEntropy(obj.weight.seqvect) ;
           if (nargin>1 && numel(reads)>0)
               obj.weight.varcov = wcoverage(obj.weight.seqvect,size(reads(1).seqvect,2),[reads.position]) ;
           end
           obj.weight.Sequence = mat2nucleo(obj.weight.seqvect) ;
       end

       function addbmu(obj)
           obj.nbmu = obj.nbmu+1 ;
       end
       
       %function split(obj)
       %end

   end
end