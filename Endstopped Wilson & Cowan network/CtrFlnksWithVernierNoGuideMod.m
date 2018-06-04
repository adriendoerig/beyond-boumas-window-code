function [Img, Data] = CtrFlnksWithVernierNoGuideMod(Cond, Height, side)

% CNTRFLNKSWITHVERNIERNOGUIDE Generate an image containing centered 
% flanks with a vernier in the middle.  Flank height is adjustable
% as is the number of flanks. modified to work with Aaron's new model
% 
% Height = 0.5(small), 1 (same size as vernier) or 2 (large)
% 
% [Img,~] = CtrFlnksWithVernierNoGuideMod(nFlank,FlankHeight,side);
% 
% figure;
% imagesc(Img); axis image; colormap gray;
% 
% Created by Aaron Clarke, EPFL, Herzog Lab, March 7, 2011 modified by A.
% Doerig, 2016

Data = [11 30 20;
    1 1 1];
clarkeFactor = 3.56;% just to check if Aaron's code works only because of the specific stim size he chose. Set to one to cancel effects.

SecPerPix = 275;                    % arc seconds per pixel
SecPerDeg = 3600;                   % seconds of arc per degree of visual angle
PixPerSec = 1/SecPerPix;            % Pixels per second of arc
PixPerMin = 60/SecPerPix;           % Pixels per minute of arc
PixPerDeg = SecPerDeg/SecPerPix;    % Pixels per degree of visual angle

Offst = 1.5*2;                % Vernier offset (pix)
ECC_FIX = 0;                % [deg] from the center
VERNIER_LEN = 2400*clarkeFactor;         % [arcsec], one segment
VERNIER_GAP = 240*clarkeFactor;          % [arcsec]
POINTER_GAP = 18000*clarkeFactor;        % [arcsec] for the pointers
ECC_VER = 0;
INT_FLANKER = 3600*clarkeFactor;         % [arcsec]
FIX_DIA = 8*clarkeFactor;                % [arcmin]
sdLine = 0.5*clarkeFactor;                           % standard error of the gaussian,
                                        % corresponds to line width
freq = 10; % spacing between flankers
%-----------

% Image dimensions
height = 210;   % 768
width = 450;    % 1024
Img = zeros(height,width);

% Set up a Cartesian coordinate system over the image
[x,y] = meshgrid(1:width,1:height);

Xcenter = 0.5*width;   % center X coordinates
Ycenter = 0.5*height;  % center Y coordinates

% Fixation point Rect (it's red in the original stimulus)
rFixation = round(Ycenter-PixPerMin*(FIX_DIA/2)):round(Ycenter+PixPerMin*(FIX_DIA/2));
cFixation = round(Xcenter-PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX):round(Xcenter+PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX);

% Vernier gap
innerVer = 0.5*PixPerSec*(VERNIER_GAP);
inner = 0;       % vernier gap (divided by 2)
outer = inner+PixPerSec*(VERNIER_LEN);     % inner+VERNIER LEN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Side = side
    
    img = zeros(height,width);
    
    % Fixation point
    %img(rFixation,cFixation) = 1;
    
    %--- Draw Vernier
    vofs = (Side*Offst);                            % Vernier offset for this image    
    xVern = round(vofs*0.5*[1 -1] + (Xcenter+PixPerDeg*(ECC_VER)));
    yVern = round(([-(1+outer) -(1+innerVer) innerVer outer] + Ycenter));
    rVern1 = yVern(1):yVern(2);
    rVern2 = yVern(3):yVern(4);
    cVern1 = xVern(1);
    cVern2 = xVern(2);
    
    img(rVern1,cVern1) = 1;
    img(rVern2,cVern2) = 1;
    
    % Save the coordinates
    XCoords = [cVern1 cVern2; cVern1 cVern2];
    YCoords = [rVern1(1) rVern2(1); rVern1(end) rVern2(end)]; 
    
    %--- Draw Pointers
%     innerP = 0.5*PixPerSec*(POINTER_GAP); % pointer gap (divided by 2)
%     outerP = innerP+PixPerSec*(VERNIER_LEN); % inner2+VERNIER LEN
%     xPtr = round(0.5+(Xcenter+PixPerDeg*(ECC_VER)));
%     yPtr = round([-outerP, -innerP, innerP, outerP]+Ycenter);
%     rPtr1 = yPtr(1):yPtr(2);
%     rPtr2 = yPtr(3):yPtr(4);
%     cPtr = xPtr;
%     
%     img(rPtr1,cPtr) = 1;
%     img(rPtr2,cPtr) = 1;
    
    %--- Choose one of the conditions 
    if Cond>1
        for iSide = [-1 +1]
        
            vofs = 0;                            % Vernier offset for this image    
            xVern = floor(vofs*0.5*[1 -1] + (Xcenter+iSide*freq+PixPerDeg*(ECC_VER)));
            yVern = floor([-outer -inner inner outer]*Height + Ycenter);
            rVern1 = [yVern(1):yVern(2)];
            rVern2 = [yVern(3):yVern(4)];
            cVern1 = xVern(1);
            cVern2 = xVern(2);

            img(rVern1,cVern1) = 1;
            img(rVern2,cVern2) = 1;
        end
        
        if Cond>2 % 8 flanks on each side
            for iFlank = 2:8
                
                for iSide = [-1 +1]
                    
                    vofs = 0;                            % Vernier offset for this image    
                    xVern = floor(vofs*0.5*[1 -1] + (Xcenter+iSide*iFlank*freq+PixPerDeg*(ECC_VER)));
                    yVern = floor([-outer -inner inner outer]*Height + Ycenter);
                    rVern1 = yVern(1):yVern(2);
                    rVern2 = yVern(3):yVern(4);
                    cVern1 = xVern(1);
                    cVern2 = xVern(2);

                    img(rVern1,cVern1) = 1;
                    img(rVern2,cVern2) = 1;
                end
            end
        end
        
        if Cond>3 % 16 flanks on each side
            for iFlank = 2:16
                
                for iSide = [-1 +1]
                    
                    vofs = 0;                            % Vernier offset for this image    
                    xVern = floor(vofs*0.5*[1 -1] + (Xcenter+iSide*iFlank*freq+PixPerDeg*(ECC_VER)));
                    yVern = floor([-outer -inner inner outer]*Height + Ycenter);
                    rVern1 = yVern(1):yVern(2);
                    rVern2 = yVern(3):yVern(4);
                    cVern1 = xVern(1);
                    cVern2 = xVern(2);

                    img(rVern1,cVern1) = 1;
                    img(rVern2,cVern2) = 1;
                end
            end
        end
    end
    Img(:,:) = img;   
    
end










% 
% 
% % Generate the Vernier
% vHeight = 15;
% dx = round(vHeight/3);
% vL = 148;
% vR = 152;
% Vernier = zeros(300,140);
% Vernier(vL,(68-vHeight+1):68) = 1;
% Vernier(vR,72:(72+vHeight-1)) = 1;
% 
% % % Generate the Vernier guides
% % vGuide = zeros(300,140);
% % vGuide(148,(68-vHeight+2)-2*vHeight:(68-vHeight+1)-vHeight) = 1;
% % vGuide(152,(72+vHeight-1)+vHeight:(72+vHeight-1)+2*vHeight-1) = 1;
% 
% % Generate the grating
% HalfFlank = round(FlankHeight/2);
% LGratingInds = vL-nFlank*dx:dx:vL-dx;
% RGratingInds = vR+dx:dx:vR+nFlank*dx;
% Grating = zeros(300,140);
% Grating(LGratingInds,[70-HalfFlank:70+HalfFlank]) = 1;
% Grating(RGratingInds,[70-HalfFlank:70+HalfFlank]) = 1;
% 
% 
% img = max(cat(3,Vernier,Grating),[],3);
% Img=img';



