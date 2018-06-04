function [Img,Data,XCoords,YCoords] = MauroOctagons(Cond,side)

% MAUROOCTAGONS Generate the stimuli from Mauro's "octagons"
% crowding experiment.
%
% Usage: [Img,Data,XCoords,YCoords]= MauroOctagons(Condition);
%
% Conditions may be any integer from 1 to 11
%
% e.g.
%
% for iCond = 1:11
%     [Img,Data,XCoords,YCoords] = MauroOctagons(iCond);
%     for iVernOffst = 1:size(Img,3)
%         figure;
%         imshow(Img(:,:,iVernOffst));
%         figure('Position',[19 619 614 296]);
%         plot(XCoords,YCoords,'k-');
%         axis equal;
%         axis([0 450 0 210]);
%         axis ij
%         set(gca,'XTick',[],'YTick',[]);        
%     end
% end
% 
% Created by Mauro Manassi (Experiment_octagons.m)
% Modified by Aaron Clarke, 14 October 2014, EPFL, LPSY Herzog Lab to only
% produce the stimuli and not run the experiment.

data = [129.5	2000	286.5;
    85.5	1989	1983.5;
    212	1978	634.5;
    87	1346.333333	765.3333333;
    419.5	2000	928.5;
    142.5	1983.5	471.3333333];
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


SPACING = 2.2;                  % [deg]  %from the center for the vernier
Ecc_Ver = [-3:-1 1:3]*SPACING + ECC_VER; % Flanker Spacing
FLANKER_LEN = 5400;             % [arcsec]

%-----------

% Fixation point Rect (it's red in the original stimulus)
rFixation = round(Ycenter-PixPerMin*(FIX_DIA/2)):round(Ycenter+PixPerMin*(FIX_DIA/2));
cFixation = round(Xcenter-PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX):round(Xcenter+PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX);


% Vernier gap
inner = 0.5*PixPerSec*VERNIER_GAP;       % vernier gap (divided by 2)
outer = inner+PixPerSec*VERNIER_LEN;     % inner+VERNIER LEN


% Baseline shape coordinates
xShape = [-1.2071 -0.5 0.5 1.2071 1.2071 0.5 -0.50 -1.2071;
    -0.5 0.5 1.2071 1.2071 0.5 -0.5 -1.2071 -1.2071];
yShape = [-0.5 -1.2071 -1.2071 -0.5 0.5 1.2071 1.2071 0.5;
    -1.2071 -1.2071 -0.5 0.5 1.2071 1.2071 0.5 -0.5];
XShape = [xShape(:)'; yShape(:)'];

% Rotations
Rot11 = 11.25 * pi/180;
Rot22 = 22.5 * pi/180;
Rot33 = 33.75 * pi/180;
Rot75 = 75 * pi/180;
Rot349 = -11.25 * pi/180;

% octagon15 = [cos(Rot11)  -sin(Rot11); sin(Rot11) cos(Rot11)]*octagon;
XShape11 = [cos(Rot11)  -sin(Rot11); sin(Rot11) cos(Rot11)]*XShape;
xShape11 = reshape(XShape11(1,:),[2 8]);
yShape11 = reshape(XShape11(2,:),[2 8]);

XShape349 = [cos(Rot349)  -sin(Rot349); sin(Rot349) cos(Rot349)]*XShape;
xShape349 = reshape(XShape349(1,:),[2 8]);
yShape349 = reshape(XShape349(2,:),[2 8]);

XShape22 = [cos(Rot22)  -sin(Rot22); sin(Rot22) cos(Rot22)]*XShape;
xShape22 = reshape(XShape22(1,:),[2 8]);
yShape22 = reshape(XShape22(2,:),[2 8]);

XShape33 = [cos(Rot33)  -sin(Rot33); sin(Rot33) cos(Rot33)]*XShape;
xShape33 = reshape(XShape33(1,:),[2 8]);
yShape33 = reshape(XShape33(2,:),[2 8]);
    
% Hexagon radius
HEX_R = FLANKER_LEN/2;    
    
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
    
    if Cond>1
        
        if Cond>=7
            % Shape around the Vernier
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*HEX_R*xShape22;
            yCoords = Ycenter+PixPerSec*HEX_R*yShape22;
            
        else
            
            % Shape around the Vernier
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*HEX_R*xShape;
            yCoords = Ycenter+PixPerSec*HEX_R*yShape;
            
        end
        
        ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
        img = max(ShapeImg,img);
        
        XCoords = cat(2,XCoords,xCoords);
        YCoords = cat(2,YCoords,yCoords);
    end
      
    if Cond==3 % 7 octagons
            
        % Make the flanking shapes
        for iEcc = 1:length(Ecc_Ver)
            
            xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*xShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end
        
    elseif Cond==4 % Normal centered and rest rotated by pi/4
            
              
        % Make the flanking shapes
        yCoordsR = Ycenter+PixPerSec*HEX_R*yShape349;
        yCoordsL = Ycenter+PixPerSec*HEX_R*yShape11;
        for iEcc = 1:length(Ecc_Ver)
            
            if Ecc_Ver(iEcc)<0
                xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*xShape11;
                yCoords = yCoordsL;
            else
                xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*xShape349;
                yCoords = yCoordsR;
            end
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end
        
        
    elseif Cond==5 % Normal centered and rest rotated by 22.5 deg
        
        % Make the flanking shapes
        yCoords = Ycenter+PixPerSec*HEX_R*yShape22;
        for iEcc = 1:length(Ecc_Ver)
            xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*xShape22;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end                  
            
            
    elseif Cond==6 % Normal centered and rest rotated by +/- 11.25 deg
            
        % Make the flanking shapes
        yCoordsL = Ycenter+PixPerSec*HEX_R*yShape349;
        yCoordsR = Ycenter+PixPerSec*HEX_R*yShape11;
        for iEcc = 1:length(Ecc_Ver)
            
            if Ecc_Ver(iEcc)<0  % Left of the Vernier
                xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*xShape349;
                yCoords = yCoordsL;
            else                % Right of the Vernier
                xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*xShape11;
                yCoords = yCoordsR;
            end
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end                      
            
            
    elseif Cond==8 % 7 octagons all rotated by 22 degrees
        
        
        % Make the flanking shapes
        yCoords = Ycenter+PixPerSec*HEX_R*yShape22;
        for iEcc = 1:length(Ecc_Ver)
            xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*xShape22;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end           
            
    elseif Cond==9 % 7 octagons, centered rotated by 90 degrees, rest converging
            
        
        % Make the flanking shapes
        yCoordsL = Ycenter+PixPerSec*HEX_R*yShape33;
        yCoordsR = Ycenter+PixPerSec*HEX_R*yShape11;
        for iEcc = 1:length(Ecc_Ver)
            
            if Ecc_Ver(iEcc)<0  % Left of the Vernier
                xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*xShape33;
                yCoords = yCoordsL;
            else                % Right of the Vernier
                xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*xShape11;
                yCoords = yCoordsR;
            end
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end      
          
            
    elseif Cond==10 % 7 octagons (centered at 90 degrees, rest normal)
            
        % Make the flanking shapes
        yCoords = Ycenter+PixPerSec*HEX_R*yShape;
        for iEcc = 1:length(Ecc_Ver)
            
            xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*xShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end
        
            
    elseif Cond==11 % 7 octagons, centered rotated by 90 degrees, rest diverging
                    
               
        % Make the flanking shapes
        yCoordsL = Ycenter+PixPerSec*HEX_R*yShape11;
        yCoordsR = Ycenter+PixPerSec*HEX_R*yShape33;
        for iEcc = 1:length(Ecc_Ver)
            
            if Ecc_Ver(iEcc)<0  % Left of the Vernier
                xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*xShape11;
                yCoords = yCoordsL;
            else                % Right of the Vernier
                xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*xShape33;
                yCoords = yCoordsR;
            end
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
        end   
                        
    end

    Img(:,:) = img;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




