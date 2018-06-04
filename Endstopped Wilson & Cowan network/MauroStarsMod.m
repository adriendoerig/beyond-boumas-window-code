function [Img,Data,XCoords,YCoords] = MauroStars(Cond,side)

% MAUROSTARS Generate the stimuli from Mauro's "stars"
% crowding experiment.
%
% Usage: [Img,Data,XCoords,YCoords] = MauroStars(Condition);
%
% Conditions may be any integer from 1 to 9
%
% e.g.
%
% for iCond = 1:9
%     [Img,Data,XCoords,YCoords] = MauroStars(iCond);
%     for iVernOffst = 1:size(Img,3)
%         figure;
%         imshow(Img(:,:,iVernOffst));        
%     end
%     figure('Position',[19 619 614 296]);
%     plot(XCoords,YCoords,'k-');
%     axis equal;
%     axis([0 450 0 210]);
%     axis ij;
%     set(gca,'XTick',[],'YTick',[]);
% end
%
% Created by Mauro Manassi (Experiment_same_stars.m)
% Modified by Aaron Clarke, 14 October 2014, EPFL, LPSY Herzog Lab to only
% produce the stimuli and not run the experiment.

data = [64.5	1991.5	535.3333333;
    79.5	1488.666667	512.6666667;
    91.5	1967.5	158.5;
    371	1385.5	734.5;
    247	1966	648];
Data = zeros(2,length(data(1,:)));
Data(1,:) = sum(data,1)/length(data(:,1));
Data(2,:) = std(data,1)/sqrt(length(data(:,1))); 

clarkeFactor = 3.56; %scale to match enStoppedModel's params.

% Stimulus conversion factors
SecPerPix = 275;                    % arc seconds per pixel
SecPerDeg = 3600;                   % seconds of arc per degree of visual angle
PixPerSec = 1/SecPerPix;            % Pixels per second of arc
PixPerMin = 60/SecPerPix;           % Pixels per minute of arc
PixPerDeg = SecPerDeg/SecPerPix;    % Pixels per degree of visual angle

Offst = 1.5*2;                        % Vernier offset (pix)
ECC_FIX = 9*clarkeFactor;                        % [deg] %from the center
VERNIER_LEN = 2400*clarkeFactor;                 % [arcsec], one segment
VERNIER_GAP = 240*clarkeFactor;                  % [arcsec]
POINTER_GAP = 18000*clarkeFactor;                % [arcsec] for the pointers
ECC_VER = 0;
FLANKER_LEN = 7200*clarkeFactor;                 % [arcsec]
IRR_R = 1.4*(FLANKER_LEN/2)*clarkeFactor;
FIX_DIA = 8*clarkeFactor;                        % [arcmin]
INT_FLANKER = 3600*clarkeFactor;                 % [arcsec]
sdLine = 0.5;                       % standard error of the gaussian,
                                    % corresponds to line width
freq = 10;                          % spacing between flankers
                                    % corresponds to line width
%-----------

% Image dimensions
height = 500;
width = 1000;
ImSz = [height width];
Img = zeros(height,width,1);

Xcenter = 0.5*width;   % center X coordinates
Ycenter = 0.5*height;  % center Y coordinates

% Fixation point Rect (it's red in the original stimulus)
rFixation = round(Ycenter-PixPerMin*(FIX_DIA/2)):round(Ycenter+PixPerMin*(FIX_DIA/2));
cFixation = round(Xcenter-PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX):round(Xcenter+PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX);


% Vernier gap
inner = 0.5*PixPerSec*VERNIER_GAP;       % vernier gap (divided by 2)
outer = inner+PixPerSec*VERNIER_LEN;     % inner+VERNIER LEN


%--- Stars coordinates

% 4-star
nPts = 4;
Amplitude = 0.3;
phi = 1.5*pi;
theta = (pi/2-phi)/nPts:pi/nPts:2*pi;
rm = 1.3225*(1+Amplitude*sin(nPts*theta+phi));
[x4,y4] = pol2cart(theta,rm);
x4 = [x4(1:end-2); x4(2:end-1)];
y4 = [y4(1:end-2); y4(2:end-1)];


% 7-star
nPts = 7;
phi = 2*pi;
theta = (pi/2-phi)/nPts:pi/nPts:2*pi;
rm = 1.3225*(1+Amplitude*sin(nPts*theta+phi));
[x7,y7] = pol2cart(theta,rm);
x7 = [x7(1:end-2); x7(2:end-1)];
y7 = [y7(1:end-2); y7(2:end-1)];

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
    
    
    %--- Choose one of the conditions from A to I
    if Cond>1 && Cond<6
        
        
        % Draw a 4-star around the Vernier
        xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x4;
        yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*y4;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
        
        if Cond==3 % 7 4-stars
            
            SPACING = 2.2*clarkeFactor;% [deg]  %from the center for the vernier
            Ecc_Ver = [-3:-1 1:3]*SPACING + ECC_VER; % Flanker Spacing                      
            
            % Make the flanking shapes
            for iEcc = 1:length(Ecc_Ver)
                
                xCoords = Xcenter+PixPerDeg*Ecc_Ver(iEcc)+PixPerSec*(INT_FLANKER)*x4;
                ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
                img = max(ShapeImg,img);
                
                XCoords = cat(2,XCoords,xCoords);
                YCoords = cat(2,YCoords,yCoords);
            end           
            
            
        elseif Cond==4 %7 4-stars touching
            
            SPACING = 2.4*clarkeFactor;% [deg]  %from the center for the vernier
            Ecc_Ver = [-3:-1 1:3]*SPACING + ECC_VER; % Flanker Spacing                      
            
            % Make the flanking shapes
            for iEcc = 1:length(Ecc_Ver)
                
                xCoords = Xcenter+PixPerDeg*Ecc_Ver(iEcc)+PixPerSec*(INT_FLANKER)*x4;
                ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
                img = max(ShapeImg,img);
                
                XCoords = cat(2,XCoords,xCoords);
                YCoords = cat(2,YCoords,yCoords);
            end  
            
            
        elseif Cond==5 %7 4-stars largest spacing (like for Mix experiment)
            
            SPACING = 2.9095*clarkeFactor; % [deg]  %from the center for the vernier
            Ecc_Ver = [-3:-1 1:3]*SPACING + ECC_VER; % Flanker Spacing                      
            
            % Make the flanking shapes
            for iEcc = 1:length(Ecc_Ver)
                
                xCoords = Xcenter+PixPerDeg*Ecc_Ver(iEcc)+PixPerSec*(INT_FLANKER)*x4;
                ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
                img = max(ShapeImg,img);
                
                XCoords = cat(2,XCoords,xCoords);
                YCoords = cat(2,YCoords,yCoords);
            end  
        end
    end
    if Cond>5
        
        % Draw a 7-star around the Vernier
        xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*INT_FLANKER*x7;
        yCoords = Ycenter+PixPerSec*(FLANKER_LEN/2)*y7;
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);                    
            
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
            
        if Cond==7 % 7 7-stars larger spacing
            
            SPACING = 2.53*clarkeFactor;% [deg]  %from the center for the vernier
            Ecc_Ver = [-3:-1 1:3]*SPACING + ECC_VER; % Flanker Spacing                      
            
            % Make the flanking shapes
            for iEcc = 1:length(Ecc_Ver)
                
                xCoords = Xcenter+PixPerDeg*Ecc_Ver(iEcc)+PixPerSec*(INT_FLANKER)*x7;
                ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
                img = max(ShapeImg,img);
                
                XCoords = cat(2,XCoords,xCoords);
                YCoords = cat(2,YCoords,yCoords);
            end  
            
        elseif Cond==8 % 7 7-stars
            
            SPACING = 2.2;% [deg]  %from the center for the vernier
            Ecc_Ver = [-3:-1 1:3]*SPACING*clarkeFactor + ECC_VER; % Flanker Spacing                      
            
            % Make the flanking shapes
            for iEcc = 1:length(Ecc_Ver)
                
                xCoords = Xcenter+PixPerDeg*Ecc_Ver(iEcc)+PixPerSec*(INT_FLANKER)*x7;
                ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
                img = max(ShapeImg,img);
                
                XCoords = cat(2,XCoords,xCoords);
                YCoords = cat(2,YCoords,yCoords);
            end 
            
            
            
        elseif Cond==9 % 7 7-stars larger spacing (same spacing as in Mix)
            
            SPACING = 2.9095*clarkeFactor;% [deg]  %from the center for the vernier
            Ecc_Ver = [-3:-1 1:3]*SPACING + ECC_VER; % Flanker Spacing                      
            
            % Make the flanking shapes
            for iEcc = 1:length(Ecc_Ver)
                
                xCoords = Xcenter+PixPerDeg*Ecc_Ver(iEcc)+PixPerSec*(INT_FLANKER)*x7;
                ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
                img = max(ShapeImg,img);
                
                XCoords = cat(2,XCoords,xCoords);
                YCoords = cat(2,YCoords,yCoords);
            end 
                        
        end
    end
    
    Img(:,:) = img;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




