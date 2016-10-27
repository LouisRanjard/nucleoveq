function [dist] = vqdist(nod,read)
% return distance between a read and a classification tree node

x = nod.weight.seqvect(1:4,read.position:(read.position+size(read.seqvect,2)-1)) ;
y = read.seqvect(1:4,:) ;

% average Euclidean distance between each column vector
dist = sum( sqrt(sum(bsxfun(@minus,x,y).^2,1)) )/size(read.seqvect,2) ;

end