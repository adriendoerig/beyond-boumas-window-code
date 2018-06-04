function [Img,Data] = MauroHexagonsMod(Cond,side)

% MAUROHEXAGONS Generate the stimuli from Mauro's "hexagons"
% crowding experiment.
%
% Usage: [Img,Data,XCoords,YCoords] = MauroHexagons(Condition);
%
% Conditions may be any integer from 1 to 11
%
% e.g.
%
% for Cond = 1:11
%     [Img,Data,XCoords,YCoords] = MauroHexagons(Cond);
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
% Created by Mauro Manassi (Experiment_hexagons.m)
% Modified by Aaron Clarke, 24 October 2014, EPFL, LPSY Herzog Lab to only
% produce the stimuli and not run the experiment.

% Human Data
data = [346.5 1791.5 1147;
    42.5 1787 543.5;
    436 1486 662;
    89 1151.5 226;
    100 1607.5 848;
    157 1747.5 1657;
    312 1321 958];
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
Img = zeros(height,width);

Xcenter = 0.5*width;   % center X coordinates
Ycenter = 0.5*height;  % center Y coordinates

SPACING = 2.2;                  % [deg]  %from the center for the vernier
Ecc_Ver = [-3:-1 1:3]*SPACING + ECC_VER; % Flanker Spacing

%-----------

% Fixation point Rect (it's red in the original stimulus)
rFixation = round(Ycenter-PixPerMin*(FIX_DIA/2)):round(Ycenter+PixPerMin*(FIX_DIA/2));
cFixation = round(Xcenter-PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX):round(Xcenter+PixPerMin*(FIX_DIA/2)-PixPerDeg*ECC_FIX);


% Vernier gap
inner = 0.5*PixPerSec*VERNIER_GAP;       % vernier gap (divided by 2)
outer = inner+PixPerSec*VERNIER_LEN;     % inner+VERNIER LEN


% Baseline shape coordinates
xShape = [-0.5 0.5 1 0.5 -0.5 -1;
    0.5 1 0.5 -0.5 -1 -0.5];
yShape = [-0.866 -0.866 0 0.866 0.866 0;
    -0.866 0 0.866 0.866 0 -0.866];

[r,c] = size(xShape);
XShape = [xShape(:)'; yShape(:)'];

% --- Rotations ---
theta = [-15 15 75 90 105]*pi/180;
nTheta = length(theta);
RotShape = zeros([r c 2 nTheta]);
for iTheta = 1:nTheta
    
    RotCoords = [cos(theta(iTheta))  -sin(theta(iTheta)); 
        sin(theta(iTheta)) cos(theta(iTheta))]*XShape;
    RotShape(:,:,1,iTheta) = reshape(RotCoords(1,:),[r c]);
    RotShape(:,:,2,iTheta) = reshape(RotCoords(2,:),[r c]);    
            
end

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
        
       
      if Cond==2 % Centered hexagon
            
            % Shape around the Vernier
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*HEX_R*xShape;
            yCoords = Ycenter+PixPerSec*HEX_R*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
              
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
      elseif Cond==3 % 7 hexagons
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*HEX_R*xShape;
            yCoords = Ycenter+PixPerSec*HEX_R*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
             
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            % Make the flanking shapes
            for iEcc = 1:length(Ecc_Ver)

                xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*xShape;
                ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
                img = max(ShapeImg,img);
                
                XCoords = cat(2,XCoords,xCoords);
                YCoords = cat(2,YCoords,yCoords);
            end             
            
      elseif Cond==4 % Normal centered and rest rotated by pi/4
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*HEX_R*xShape;
            yCoords = Ycenter+PixPerSec*HEX_R*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);   
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            % Make the flanking shapes
            
            for iEcc = 1:length(Ecc_Ver)
                if Ecc_Ver(iEcc)>0
                    xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*RotShape(:,:,1,1);
                    yCoords = Ycenter+PixPerSec*HEX_R*RotShape(:,:,2,1);
                else
                    xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*RotShape(:,:,1,2);
                    yCoords = Ycenter+PixPerSec*HEX_R*RotShape(:,:,2,2);
                end
                ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
                img = max(ShapeImg,img);
                
                XCoords = cat(2,XCoords,xCoords);
                YCoords = cat(2,YCoords,yCoords);
            end    
                  
          
      elseif Cond==5 % Normal centered and rest rotated by pi/2
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*HEX_R*xShape;
            yCoords = Ycenter+PixPerSec*HEX_R*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);                         
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            % Make the flanking shapes
            yCoords = Ycenter+PixPerSec*HEX_R*RotShape(:,:,2,4);
            for iEcc = 1:length(Ecc_Ver)

                xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*RotShape(:,:,1,4);
                ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
                img = max(ShapeImg,img);
                
                XCoords = cat(2,XCoords,xCoords);
                YCoords = cat(2,YCoords,yCoords);
            end 
            
      elseif Cond==6 % Normal centered and rest rotated by -pi/4
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*HEX_R*xShape;
            yCoords = Ycenter+PixPerSec*HEX_R*yShape;
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);                         
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            % Make the flanking shapes            
            for iEcc = 1:length(Ecc_Ver)
                if Ecc_Ver(iEcc)>0
                    xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*RotShape(:,:,1,2);
                    yCoords = Ycenter+PixPerSec*HEX_R*RotShape(:,:,2,2);
                else
                    xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*RotShape(:,:,1,1);
                    yCoords = Ycenter+PixPerSec*HEX_R*RotShape(:,:,2,1);
                end
                
                ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
                img = max(ShapeImg,img);
                
                XCoords = cat(2,XCoords,xCoords);
                YCoords = cat(2,YCoords,yCoords);
            end
            
      elseif Cond==7 % Centered hexagon 90 degrees rotated
              
              xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*HEX_R*RotShape(:,:,1,4);
              yCoords = Ycenter+PixPerSec*HEX_R*RotShape(:,:,2,4);
              ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
              img = max(ShapeImg,img);
              
              XCoords = cat(2,XCoords,xCoords);
              YCoords = cat(2,YCoords,yCoords);
               
      elseif Cond==8 % 7 hexagons all rotated by 90 degrees            
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*HEX_R*RotShape(:,:,1,4);
            yCoords = Ycenter+PixPerSec*HEX_R*RotShape(:,:,2,4);
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
             
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            % Make the flanking shapes            
            for iEcc = 1:length(Ecc_Ver)

                xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*RotShape(:,:,1,4);
                ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
                img = max(ShapeImg,img);
                
                XCoords = cat(2,XCoords,xCoords);
                YCoords = cat(2,YCoords,yCoords);
            end                                     
                    
      elseif Cond==9 % 7 hexagons, centered rotated by 90 degrees, rest converging           
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*HEX_R*RotShape(:,:,1,4);
            yCoords = Ycenter+PixPerSec*HEX_R*RotShape(:,:,2,4);
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            % Make the flanking shapes            
            for iEcc = 1:length(Ecc_Ver)
                if Ecc_Ver(iEcc)>0
                    xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*RotShape(:,:,1,3);
                    yCoords = Ycenter+PixPerSec*HEX_R*RotShape(:,:,2,3);
                else
                    xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*RotShape(:,:,1,5);
                    yCoords = Ycenter+PixPerSec*HEX_R*RotShape(:,:,2,5);
                end
                ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
                img = max(ShapeImg,img);
                
                XCoords = cat(2,XCoords,xCoords);
                YCoords = cat(2,YCoords,yCoords);
            end                
                    
      elseif Cond==10 % 7 hexagons (centered at 90 degrees, rest normal)            
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*HEX_R*RotShape(:,:,1,4);
            yCoords = Ycenter+PixPerSec*HEX_R*RotShape(:,:,2,4);
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
            
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            % Make the flanking shapes 
            yCoords = Ycenter+PixPerSec*HEX_R*yShape;
            for iEcc = 1:length(Ecc_Ver)

                xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*xShape;
                ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
                img = max(ShapeImg,img);
                
                XCoords = cat(2,XCoords,xCoords);
                YCoords = cat(2,YCoords,yCoords);
            end               
            
      elseif Cond==11 % 7 hexagons, centered rotated by 90 degrees, rest diverging                        
            
            xCoords = Xcenter+PixPerDeg*ECC_VER+PixPerSec*HEX_R*RotShape(:,:,1,4);
            yCoords = Ycenter+PixPerSec*HEX_R*RotShape(:,:,2,4);
            ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
            img = max(ShapeImg,img);
                        
            XCoords = cat(2,XCoords,xCoords);
            YCoords = cat(2,YCoords,yCoords);
            
            % Make the flanking shapes            
            for iEcc = 1:length(Ecc_Ver)
                if Ecc_Ver(iEcc)>0
                    xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*RotShape(:,:,1,5);
                    yCoords = Ycenter+PixPerSec*HEX_R*RotShape(:,:,2,5);
                else
                    xCoords = Xcenter + PixPerDeg*Ecc_Ver(iEcc) + PixPerSec*HEX_R*RotShape(:,:,1,3);
                    yCoords = Ycenter+PixPerSec*HEX_R*RotShape(:,:,2,3);
                end
                ShapeImg = AntiAlsLine(xCoords,yCoords,ImSz,sdLine);
                img = max(ShapeImg,img);
                
                XCoords = cat(2,XCoords,xCoords);
                YCoords = cat(2,YCoords,yCoords);
            end 
            
      end
      Img(:,:) = img;
end

     
        
        
        

       
    
   
    

    