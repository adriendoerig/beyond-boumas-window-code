function Img = AntiAlsLine(xCoords,yCoords,Sz,S)

% ANTIALSLINE Generate an image containing lines using Gaussian anti-aliasing.
% 
% Usage: Img = AntiAlsLine(xCoords,yCoords,Sz,S);
% 
% e.g.
% 
% % Standard deviation on Gaussian around line
% S = 1;
% 
% % Image size [rows columns]
% Sz = [600 800];
% 
% % Line end coordinates
% xCoords = [100 200 300;
%     150 200 300];
% yCoords = [1 50 25;
%     50 100 75];
% 
% Img = AntiAlsLine(xCoords,yCoords,Sz,S);
% 
% figure;
% imagesc(Img);
% axis image; 
% colormap gray;
% 
% Created by Aaron Clarke, EPFL, Herzog Lab, July 25, 2011

% Number of lines
nLine = size(xCoords,2);

% Cartesian coordinate system over the image
x = 1:Sz(2);
y = 1:Sz(1);

% Imaging lines
[X,Y] = meshgrid(x,y);

theta = atan2(yCoords(2,:)-yCoords(1,:),xCoords(2,:)-xCoords(1,:));
w1 = sin(theta);
w2 = -cos(theta);
w0 = -w1.*xCoords(1,:) - w2.*yCoords(1,:);

Dists = repmat(reshape(w1,[1 1 nLine]),[Sz 1]).*repmat(X,[1 1 nLine]) + ...
    repmat(reshape(w2,[1 1 nLine]),[Sz 1]).*repmat(Y,[1 1 nLine]) + ...
    repmat(reshape(w0,[1 1 nLine]),[Sz 1]);

% Exponential function of line distance
TheExp = exp(-Dists.^2/(2*S^2));

% Calculate the orthoganal distances from the line end-points
OrthTheta = theta+pi/2;
w1Orth = sin(OrthTheta);
w2Orth = -cos(OrthTheta);
for iCoord = 1:2
    w0orth = -w1Orth.*xCoords(iCoord,:) - w2Orth.*yCoords(iCoord,:);
    OrthDist = repmat(reshape(w1Orth,[1 1 nLine]),[Sz 1]).*repmat(X,[1 1 nLine]) + ...
        repmat(reshape(w2Orth,[1 1 nLine]),[Sz 1]).*repmat(Y,[1 1 nLine]) + ...
        repmat(reshape(w0orth,[1 1 nLine]),[Sz 1]);

    for iLine = 1:nLine
        OtherInd = 1-iCoord+2;
        TheSign = sign(w1Orth(iLine)*xCoords(OtherInd,iLine) + w2Orth(iLine)*yCoords(OtherInd,iLine) + w0orth(iLine));
        OrthDist(:,:,iLine) = TheSign*OrthDist(:,:,iLine);
    end
    TheExp(OrthDist<0) = 0;    
end


% Combine the lines to make the final image
% Img = round(max(TheExp,[],3)*254+1);
Img = max(TheExp,[],3);

