function [Img,Data] = MauroCirclesMod(Cond,side)

% Make the Vernier with crowding circles images from Mauro's experiments.
%
% Usage: [Img,Data,XCoords,YCoords] = MauroCircles(Cond);
%
% Cond may be any integer from 1 to 6.
%
% e.g.
%
% for iCond = 1:6
%     [Img,Data,XCoords,YCoords] = MauroCircles(iCond);
%     for iVernOffst = 1:size(Img,3)
%         figure;
%         imshow(Img(:,:,iVernOffst));
%         figure('Position',[19 619 614 296]);
%         plot(XCoords,YCoords,'k-');
%         axis equal;
%         axis([0 450 0 210]);
%         axis ij;
%         set(gca,'XTick',[],'YTick',[]);
%     end
% end
%
% Created by Mauro Manassi, EPFL, LPSY, Herzog Lab (Experiment_circle.m)
% Modified by Aaron Clarke Oct. 7, 2014 to return only the image and not
% run the entire expeirment.

% Circles with different spacing

data = [281.3333333	825	455.5;
    158	1079	557.6666667;
    75.5	1397	501;
    133.5	1998	425.5];
Data = zeros(2,3);
Data(1,:) = sum(data,1)/length(data(:,1));
Data(2,:) = std(data,1)/sqrt(length(data(:,1)));

clarkeFactor = 3.56; % set size to match aaron's params in endStoppedModel
% Stimulus conversion factors

SecPerPix = 275;                    % arc seconds per pixel
SecPerDeg = 3600;                   % seconds of arc per degree of visual angle
PixPerSec = 1/SecPerPix;            % Pixels per second of arc
PixPerMin = 60/SecPerPix;           % Pixels per minute of arc
PixPerDeg = SecPerDeg/SecPerPix;    % Pixels per degree of visual angle

Offst = 1.5*2;                % Vernier offset (pix)
ECC_FIX = 9*clarkeFactor;                % [deg] from the center
VERNIER_LEN = 2400*clarkeFactor;         % [arcsec], one segment
VERNIER_GAP = 240*clarkeFactor;          % [arcsec]
POINTER_GAP = 18000*clarkeFactor;        % [arcsec] for the pointers
ECC_VER = 0*clarkeFactor;
INT_FLANKER = 3600*clarkeFactor;         % [arcsec]
FIX_DIA = 8*clarkeFactor;                % [arcmin]
CircRad = PixPerSec*INT_FLANKER*3/2.3;  % Circle Radius
sdLine = 0.5;                           % standard error of the gaussian,
                                        % corresponds to line width
FlnkDists = [0 0 PixPerSec*(7200+1000) 2*CircRad/clarkeFactor PixPerSec*(7200+3000) PixPerSec*(7200+4000)].*clarkeFactor;                                         
nEllipsePoints = 50;
%-----------

% Image dimensions
height = 500;   % 768
width = 1000;    % 1024
Img = zeros(height,width);

% Set up a Cartesian coordinate system over the image
[x,y] = meshgrid(1:width,1:height);

Xcenter = 0.5*width;   % center X coordinates
Ycenter = 0.5*height;  % center Y coordinates

% Fixation point Rect (it's red in the original stimulus)
rFixation = round(Ycenter-PixPerMin*(FIX_DIA/2)):round(Ycenter+PixPerMin*(FIX_DIA/2));
cFixation = round(Xcenter-PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX):round(Xcenter+PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX);

% Vernier gap
inner = 0.5*PixPerSec*(VERNIER_GAP);       % vernier gap (divided by 2)
outer = inner+PixPerSec*(VERNIER_LEN);     % inner+VERNIER LEN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Side = side
    
    img = zeros(height,width);
    
    % Fixation point
    img(rFixation,cFixation) = 1;
    
    %--- Draw Vernier
    vofs = (Side*Offst);                            % Vernier offset for this image    
    xVern = round(vofs*0.5*[1 -1] + (Xcenter+PixPerDeg*(ECC_VER)));
    yVern = round([-outer -inner inner outer] + Ycenter);
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
    innerP = 0.5*PixPerSec*(POINTER_GAP); % pointer gap (divided by 2)
    outerP = innerP+PixPerSec*(VERNIER_LEN); % inner2+VERNIER LEN
    xPtr = round(0.5+(Xcenter+PixPerDeg*(ECC_VER)));
    yPtr = round([-outerP, -innerP, innerP, outerP]+Ycenter);
    rPtr1 = yPtr(1):yPtr(2);
    rPtr2 = yPtr(3):yPtr(4);
    cPtr = xPtr;
    
    img(rPtr1,cPtr) = 1;
    img(rPtr2,cPtr) = 1;
    
    %--- Choose one of the conditions from A to c
    if Cond>1
        
        % Make the circle around the Vernier
        
        % Calculate the distance of each point from the circle
        rad = sqrt((x-Xcenter).^2 + (y-Ycenter).^2);
        
        % Make the circle a Gaussian function of distance from the chosen radius
        d = exp(-(rad-CircRad).^2/(2*sdLine^2)); %multiply by 255 to get the psychtoolbox colors
        
        img = max(img,d);
        
        [cx,cy] = ellipse(nEllipsePoints,CircRad,CircRad,Xcenter,Ycenter,0); 
        XCoords = cat(2,XCoords,[cx(1:end-1); cx(2:end)]);
        YCoords = cat(2,YCoords,[cy(1:end-1); cy(2:end)]);
        
        if Cond>2 % 7 circles
                        
            % Make the flanking circles
            for iCirc = 1:3
                
                for iSide = [-1 +1]
                    
                    % Radius from circle center
                    rad = sqrt((x-Xcenter+iSide*iCirc*FlnkDists(Cond)).^2 + (y-Ycenter).^2);
                    
                    % Gaussian function of radius
                    d = exp(-(rad-CircRad).^2/(2*sdLine^2));
                    
                    img = max(img,d);
                    
                    [cx,cy] = ellipse(nEllipsePoints,CircRad,CircRad,...
                        Xcenter+iSide*iCirc*FlnkDists(Cond),Ycenter,0); 
                    XCoords = cat(2,XCoords,[cx(1:end-1); cx(2:end)]);
                    YCoords = cat(2,YCoords,[cy(1:end-1); cy(2:end)]);
                end
            end
        end
    end
    Img(:,:) = img;   
    
end





