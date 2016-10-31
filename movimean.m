function [ y ] = movimean(x,windowSize)
% perform moving mean of input signal x around window size
%

if nargin<2
    windowSize = 3 ;
end

b = (1/windowSize)*ones(1,windowSize) ;
a = 1 ;
y = filter(b,a,x) ;

end