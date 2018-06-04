function [y] = iseven(x)

% ISEVEN Returns true if the argument is even and false otherwise.
%
% Usage: [y] = iseven(x);
%
% e.g.
%
% iseven(3)
% iseven(4)
%
% Created by Aaron Clarke, NDSU, for Mark McCourt, May 29, 2007
   

y = ~mod(x,2);