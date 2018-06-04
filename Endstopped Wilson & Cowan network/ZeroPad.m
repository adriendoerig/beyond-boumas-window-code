function [y] = ZeroPad(x,Sz)

% ZEROPAD Pad matrix with zeros so that it reaches the desired size.
% 
% e.g.
% 
% x = ones(2);
% Sz = [4 4];
% [y] = ZeroPad(x,Sz)
% 
% Created by Aaron Clarke, EPFL, LPSY, Herzog Lab, December 12, 2011

% Get the size of x
xSz = size(x);
if length(xSz)<length(Sz)
    xSz = [xSz, zeros(1,length(Sz)-length(xSz))];
end

% Initialize y
y = x;

for iDim = 1:length(Sz)
    
    if xSz(iDim)<Sz(iDim)
        SzDiff = Sz(iDim)-xSz(iDim);
        if iseven(SzDiff)
            Hlf = SzDiff/2;
            TheseZeros = xSz;
            TheseZeros(iDim) = Hlf;
            zInd = find(TheseZeros==0);
            if ~isempty(zInd)
                TheseZeros = TheseZeros(1:zInd(1)-1);                            
            end
            if xSz(iDim)
                y = cat(iDim,zeros(TheseZeros),y,zeros(TheseZeros));
            else
                TheseZeros2 = TheseZeros;
                TheseZeros2(iDim) = TheseZeros(iDim)-1;
                y = cat(iDim,zeros(TheseZeros),y,zeros(TheseZeros2));
            end
            
        else
            Hlf = (SzDiff-1)/2;
            TheseZeros = xSz;
            zInd = find(TheseZeros==0);
            if ~isempty(zInd)
                TheseZeros = TheseZeros(1:zInd(1)-1);
            end
            TheseZeros1 = TheseZeros;
            TheseZeros2 = TheseZeros;
            TheseZeros1(iDim) = Hlf;
            TheseZeros2(iDim) = Hlf+1;
            if xSz(iDim)
                y = cat(iDim,zeros(TheseZeros1),y,zeros(TheseZeros2));
            else                
                TheseZeros2(iDim) = TheseZeros2(iDim)-1;
                y = cat(iDim,zeros(TheseZeros),y,zeros(TheseZeros2));
            end
        end
    end
    
    % Re-calculate xSz
    xSz = size(y);
    if length(xSz)<length(Sz)
        xSz = [xSz, zeros(1,length(Sz)-length(xSz))];
    end        
        
end
        
    