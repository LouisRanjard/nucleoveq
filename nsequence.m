classdef nsequence
% write a description of the class here.

   properties
        Header ;
        Sequence ;
        Quality ;
        Seqvect ;
        AliPosition ;
   end

   methods

       function obj = nsequence()
       % class constructor
           if(nargin > 0)
                obj.Header = '' ;
                obj.Sequence = 'ACTG' ;
                obj.Quality = '' ;
                obj.Seqvect = [ 1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1 ] ;
                obj.AliPosition = 0 ;
           end
       end

       function align(obj,read)
       end

   end
end