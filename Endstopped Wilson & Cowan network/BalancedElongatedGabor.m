function [G] = BalancedElongatedGabor(lambda,theta,phi)

% BALANCEDELONGATEDGABOR Produces a balanced Gabor function (i.e. a Gabor
% who's mean amplitude is 0).
% 
% Usage: [G] = BalancedElongatedGabor(Period,Orientation,Phase);
% 
% e.g.
% 
% lambda = 32;
% theta = pi/2;
% phi = 0; %pi/2;
% 
% [G] = BalancedElongatedGabor(lambda,theta,phi);
% 
% figure;
% imagesc(G); axis image; colormap gray;
% 
% Created by Aaron Clarke, NDSU, Visual Dynamics lab, February 26, 2007

% Standard deviation on the Gaussian envelope
s = lambda*2;
S = s*[2 0; 0 1];

% Rotation matrix
R = [cos(theta) sin(theta); 
    -sin(theta) cos(theta)];

% TheMax = ceil(3*s);
TheMax = ceil(lambda);
[x,y] = meshgrid(-TheMax:TheMax);   % Cartesian coordinate system

% Rotate the coordinate system
X = (R*[x(:)'; y(:)'])';

% Gaussian envelope
P = reshape(NormPDFN(X,[0 0],S),size(x));

G = P .* ...
    ( cos( 2*pi/lambda * (x*cos(theta) + y*sin(theta)) - phi) - ...
    cos(phi)*exp(-2*pi^2*(s^2)/(lambda^2)) );


% figure;
% imagesc(G); axis image; colormap gray;