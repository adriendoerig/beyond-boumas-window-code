function [Vernier,img] = JitFlnksWithVernier(nFlank)

% JITRFLNKSWITHVERNIER Generate an image containing vertically jittered 
% flanks with a vernier in the middle.  Flank height is adjustable
% as is the number of flanks.
% 
% Usage: [Vernier,img] = JitFlnksWithVernier(nFlank);
% 
% e.g.
% 
% nFlank = 8;       % Options are: 2, 4, 8 or 16
% 
% [Vernier,img] = JitFlnksWithVernier(nFlank);
% 
% figure;
% subplot(2,1,1);
% imagesc(Vernier'); axis image; colormap gray;
% subplot(2,1,2);
% imagesc(img'); axis image; colormap gray;
% 
% Created by Aaron Clarke, EPFL, Herzog Lab, March 7, 2011
% Modified by Aaron Clarke, August 11 to exactly reproduce Mauro's 
% jittered stimuli


% Vernier/Image parameters
ImSz = [450 210];% [300 140]; %                  % Image size
vHeight = 15;                       % Vernier height         
Mid = round(ImSz/2);                % Middle image coordinates
Offst = 1.5;
mid = (ImSz(1)+1)/2;
vL = round(mid-Offst);              % Vernier left position
vR = round(mid+Offst);              % Vernier right position
vGap = 2;                           % Gap between Vernier segments and midline
hGap = 5;                           % Gap between flankers
fGap = 7;

% Initialize the image
img = zeros(ImSz);

% Generate the Vernier
img(vL,(Mid(2)-vGap-vHeight+1):Mid(2)-vGap) = 1;
img(vR,Mid(2)+vGap:(Mid(2)+vGap+vHeight-1)) = 1;
Vernier = img;

Const = 0:nFlank-1; Const = [-fliplr(Const) Const]; 
cSign = ones(1,nFlank); cSign = [-cSign cSign];
Gaps = Const*hGap + cSign*fGap;
 
JittMid = []; 
if nFlank == 0
    JittMid = [];    
elseif nFlank==2
    JittMid = [-3 -1 1 -3];      
elseif nFlank==4
    JittMid = [-5 7 -3 1 1 -3 7 -5];              
elseif nFlank==8
    JittMid = [-7 -1 -6 8 -5 7 -3 1 1 -3 7 -5 8 -6 -1 -7];
elseif nFlank==16    
    JittMid = [-7 -1 -6 8 -5 7 -3 1 1 -3 7 -5 8 -6 -1 -7];
    JittMid = [3 5 5 -3 3 -8 4 -11 JittMid -11 4 -8 3 -3 5 5 3];    
end

if ~isempty(JittMid)
    for iLine = 1:length(JittMid)
        img(Mid(1)+Gaps(iLine),round(Mid(2)+JittMid(iLine)-(vHeight+1)/2):round(Mid(2)+JittMid(iLine)+(vHeight+1)/2)) = 1;
    end
end

