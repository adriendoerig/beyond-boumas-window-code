function [Img,Data] = MauroIrregular1Mod(Cond,side)

% MAUROIRREGULAR1 Generate the stimuli from Mauro's first "irregular shapes"
% crowding experiment.
%
% Usage: [Img,Data,XCoords,YCoords] = MauroIrregular1(Condition);
%
% e.g.
%
% for iCond = 1:10
%     [Img,Data,XCoords,YCoords] = MauroIrregular1(iCond);
%     for iVernOffst = 1:size(Img,3)
%         figure;
%         imshow(Img(:,:,iVernOffst));
%         figure('Position',[19 619 614 296]);
%         plot(XCoords,YCoords,'k-');
%         axis equal;
%         axis([0 450 0 210]);
%         set(gca,'XTick',[],'YTick',[]);
% 
%     end
% end
%
% Created by Mauro Manassi (Experiment_irregular.m)
% Modified by Aaron Clarke, 8 Sept 2014, EPFL, LPSY Herzog Lab to only
% produce the stimuli and not run the experiment.


% Human Data
data = [87	1976	1266;
    338.5	1497.333333	1272.5;
    151	1073.333333	505;
    141.5	1607.333333	795.3333333;
    98.5	1033.666667	471;
    180	1405.666667	851.25];
Data = zeros(2,3);
Data(1,:) = sum(data,1)/length(data(:,1));
Data(2,:) = std(data,1)/sqrt(length(data(:,1)));
    
 
% Stimulus conversion factors
SecPerPix = 275;                    % arc seconds per pixel
SecPerDeg = 3600;                   % seconds of arc per degree of visual angle
PixPerSec = 1/SecPerPix;            % Pixels per second of arc
PixPerMin = 60/SecPerPix;           % Pixels per minute of arc
PixPerDeg = SecPerDeg/SecPerPix;    % Pixels per degree of visual angle

Offst = 1.5;                        % Vernier offset (pix)
ECC_FIX = 9;                        % [deg] %from the center
VERNIER_LEN = 2400;                 % [arcsec], one segment
VERNIER_GAP = 240;                  % [arcsec]
POINTER_GAP = 18000;                % [arcsec] for the pointers
ECC_VER = 0;
FLANKER_LEN = 7200;                 % [arcsec]
IRR_R = 0.65*(FLANKER_LEN/2);
FIX_DIA = 8;                        % [arcmin]
sdLine = 0.5;                       % standard error of the gaussian,
                                    % corresponds to line width
%-----------

% Image dimensions
height = 210;
width = 450;
ImSz = [height width];
Img = zeros(height,width,1);

Xcenter = 0.5*width;   % center X coordinates
Ycenter = 0.5*height;  % center Y coordinates

% Fixation point Rect (it's red in the original stimulus)
rFixation = round(Ycenter-PixPerMin*(FIX_DIA/2)):round(Ycenter+PixPerMin*(FIX_DIA/2));
cFixation = round(Xcenter-PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX):round(Xcenter+PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX);

% Baseline shape coordinates
xIrreg = [-0.8 0.6 1.2 0.2 -1 -1.5;
    0.6 1.2 0.2 -1 -1.5 -0.8];
yIrreg = [-1 -1.7 1 1.7 0.8 1.5;
    -1.7  1 1.7 0.8 1.5 -1];

% Flipped x-coordinates
xIrregF = [0.5 -0.9 -1.5 -0.5 0.7 1.2
    -0.9 -1.5 -0.5 0.7 1.2 0.5];

% Inverted figure
IrregU = [cos(pi) -sin(pi); sin(pi) cos(pi)]*[xIrreg(:)'; yIrreg(:)'];
xIrregU = reshape(IrregU(1,:),[2 6]);
yIrregU = reshape(IrregU(2,:),[2 6]);

% Some conditions only vary the spacing between flankers
SpacingCond = [3 4 5 10];
Spacings = [2 2.2 2.53 2.2];

% Vernier gap
inner = 0.5*PixPerSec*(VERNIER_GAP);       % vernier gap (divided by 2)
outer = inner+PixPerSec*(VERNIER_LEN);     % inner+VERNIER LEN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for Side = side
    
    img = zeros(height,width);
    
    % Fixation point
    img(rFixation,cFixation) = 1;
    
    %--- Draw Vernier
    vofs = (Side*Offst);                       % Vernier offset for this image        
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
    
    
    %case 'A' %Vernier alone
    if Cond>1
        
        % Shape around the Vernier
        xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*IRR_R*xIrreg;
        yCoords = Ycenter+PixPerSec*IRR_R*yIrreg;
        
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
    end    
        
    % Conditions that vary the spacing [3 4 5 10]
    [~,ind] = ismember(Cond,SpacingCond);    
    if ind
        
        Ecc_Ver = [-3:-1 1:3]*Spacings(ind) + ECC_VER; % Flanker Spacing  
        
        % Make the flanking shapes
        for iEcc = 1:length(Ecc_Ver)
            
            xCoords = Xcenter+0.85*PixPerDeg*Ecc_Ver(iEcc)+PixPerSec*(IRR_R)*xIrreg;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end
    end
    
    % Other conditions that vary shape orientations     
    if Cond==6 % 'F' - 7 irregular forms (100%) Mirror alternating
        
        SPACING = 2.53; % [deg]  %from the center for the vernier
        Ecc_Ver = [-3:-1 1:3]*SPACING + ECC_VER; % Flanker Spacing                
        
        % Make the flanking shapes
        for iEcc = 1:length(Ecc_Ver)
            
            if abs(Ecc_Ver(iEcc))==abs(Ecc_Ver(2))
                xCoords = Xcenter+0.85*PixPerDeg*Ecc_Ver(iEcc)+PixPerSec*(IRR_R)*xIrreg;
            else
                xCoords = Xcenter+0.85*PixPerDeg*Ecc_Ver(iEcc)+PixPerSec*(IRR_R)*xIrregF;
            end
            
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end
                
        
    elseif Cond==7 % 'G' - 7 irregular forms (100%) Mirror
        
        SPACING = 2.53; % [deg]  %from the center for the vernier
        Ecc_Ver = [-3:-1 1:3]*SPACING + ECC_VER; % Flanker Spacing                
        
        % Make the flanking shapes
        for iEcc = 1:length(Ecc_Ver)
            
            if iEcc<4
                xCoords = Xcenter+0.85*PixPerDeg*Ecc_Ver(iEcc)+PixPerSec*(IRR_R)*xIrreg;
            else
                xCoords = Xcenter+0.85*PixPerDeg*Ecc_Ver(iEcc)+PixPerSec*(IRR_R)*xIrregF;
            end
            
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end
        
        
    elseif Cond==8 % 'H' - rotated 180 degrees
        
        SPACING = 2.2; % [deg]  %from the center for the vernier
        Ecc_Ver = [-3:-1 1:3]*SPACING + ECC_VER; % Flanker Spacing                
        
        % Make the flanking shapes
        yCoords = Ycenter+PixPerSec*IRR_R*yIrregU;
        for iEcc = 1:length(Ecc_Ver)
            
            xCoords = Xcenter+0.85*PixPerDeg*Ecc_Ver(iEcc)+PixPerSec*(IRR_R)*xIrregU;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end
        
        
    elseif Cond==9 % 'I' - Mirrored flankers
        
        SPACING = 2.53; % [deg]  %from the center for the vernier
        Ecc_Ver = [-3:-1 1:3]*SPACING + ECC_VER; % Flanker Spacing
                
        % Make the flanking shapes
        for iEcc = 1:length(Ecc_Ver)
            
            xCoords = Xcenter+0.85*PixPerDeg*Ecc_Ver(iEcc)+PixPerSec*(IRR_R)*xIrregF;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end
        
        
    elseif Cond==10 % 'J' - 7 irregular forms (Same spacing as in mix)
        
        SPACING = 2.2; % [deg]  %from the center for the vernier
        Ecc_Ver = [-3:-1 1:3]*SPACING + ECC_VER; % Flanker Spacing                
        
        % Make the flanking shapes
        for iEcc = 1:length(Ecc_Ver)
            
            xCoords = Xcenter+0.85*PixPerDeg*Ecc_Ver(iEcc)+PixPerSec*(IRR_R)*xIrreg;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end    
    end
    
    Img(:,:) = img;    
end







