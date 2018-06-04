function [img] = GestaltStim(StimType, side)

% GESTALTSTIM Make the Gestalt stimuli of Sayim, Westheimer and Herzog (2010).
% side: 1 for vernier to the left -1 for vernier to the right.
% 
% Usage: [img] = GestaltStim(StimType);
% 
% StimType can be either:
% 'Vernier'
% 'Lines'
% 'Rectangles'
% 'Cuboids'
% 'Scrambled Cuboids'
% '|x|'
% 'Rectangles |x|'
% 
% e.g.
% 
% [img] = GestaltStim('Cuboids');
% figure;
% imagesc(img);
% axis image;
% colormap gray;
% 
% Created by Aaron Clarke, EPFL, Herzog Lab, November 17, 2011

% Variables
clarkeFactor = 2;
ImSz = [210 450]; %[140 300];       % Image size. Use [210 450] for FeefForwardFeedbackScaledStim
vHeight = 15*clarkeFactor;                       % Vernier height
Mid = round(ImSz/2);                % Middle image coordinates
mid = (ImSz+1)/2;
Offst = 1.5;
vL = round(mid(2)-Offst);           % Vernier left position
vR = round(mid(2)+Offst);           % Vernier right position
vGap = 2;                           % Gap between Vernier segments and midline
FlnkSep = 10;                        % Separation between flankers and Vernier
FlnkHeight = 30*clarkeFactor;                    % Flanker height
RecWidth = 65*clarkeFactor;
S = 0.5;                            % Standard deviation on anti-aliased pixels
DeVectorized = false;

% Constants
MidFlnkHeight = round(FlnkHeight/2);

% Initialize the image
img = zeros(ImSz);

% Generate the Vernier
img((Mid(1)-vGap+1-vHeight+1):Mid(1)-vGap,round(mid(2)+side*Offst)) = 1;     % Top part
img(Mid(1)+vGap:(Mid(1)+vGap-1+vHeight-1),round(mid(2)-side*Offst)) = 1;     % Bottom part

if DeVectorized
    figure('Color','k');
    plot([vL vL],[(Mid(1)-vGap+1-vHeight+1) Mid(1)-vGap],'w-');
    hold on;
    plot([vR vR],[Mid(1)+vGap (Mid(1)+vGap-1+vHeight-1)],'w-');
    set(gca,'Color','k','XColor','w','YColor','w','XTick',[],'YTick',[]);
    axis([1 ImSz(2) 1 ImSz(1)]);
    axis equal;
    axis ij;
end


if ~strcmp(lower(StimType),'vernier')

    % Make flanks
    img(Mid(1)-MidFlnkHeight:Mid(1)+MidFlnkHeight,vL-FlnkSep) = 1;
    img(Mid(1)-MidFlnkHeight:Mid(1)+MidFlnkHeight,vR+FlnkSep) = 1;
    
    if DeVectorized
        plot([vL-FlnkSep vL-FlnkSep],[Mid(1)-MidFlnkHeight Mid(1)+MidFlnkHeight],'w-');
        plot([vR+FlnkSep vR+FlnkSep],[Mid(1)-MidFlnkHeight Mid(1)+MidFlnkHeight],'w-');
    end

end

switch lower(StimType)
    case 'vernier'
    case 'lines'
    case 'rectangles'
        
        EndL = vL-FlnkSep-RecWidth;
        EndR = vR+FlnkSep+RecWidth;
        
        % Outer vertical rectangle segments
        img(Mid(1)-MidFlnkHeight:Mid(1)+MidFlnkHeight,EndL) = 1;
        img(Mid(1)-MidFlnkHeight:Mid(1)+MidFlnkHeight,EndR) = 1;
        
        % Top connecting rectangle segments
        img(Mid(1)+MidFlnkHeight,EndL:vL-FlnkSep) = 1;
        img(Mid(1)+MidFlnkHeight,vR+FlnkSep:EndR) = 1;
        
        % Bottom connecting rectangle segments
        img(Mid(1)-MidFlnkHeight,EndL:vL-FlnkSep) = 1;
        img(Mid(1)-MidFlnkHeight,vR+FlnkSep:EndR) = 1;
        
        if DeVectorized
            plot([EndL EndL],[Mid(1)-MidFlnkHeight Mid(1)+MidFlnkHeight],'w-');
            plot([EndR EndR],[Mid(1)-MidFlnkHeight Mid(1)+MidFlnkHeight],'w-');
            
            plot([EndL vL-FlnkSep],[Mid(1)+MidFlnkHeight Mid(1)+MidFlnkHeight],'w-');
            plot([vR+FlnkSep EndR],[Mid(1)+MidFlnkHeight Mid(1)+MidFlnkHeight],'w-');
            
            plot([EndL vL-FlnkSep],[Mid(1)-MidFlnkHeight Mid(1)-MidFlnkHeight],'w-');
            plot([vR+FlnkSep EndR],[Mid(1)-MidFlnkHeight Mid(1)-MidFlnkHeight],'w-');
        end
        
    case 'cuboids'
        
        % Rectangle part ---------------------------------
        EndL = vL-FlnkSep-RecWidth;
        EndR = vR+FlnkSep+RecWidth;
        EndT = Mid(1)-MidFlnkHeight*clarkeFactor/2;
        EndB = Mid(1)+MidFlnkHeight*clarkeFactor/2;
        
        % Outer vertical rectangle segments
        img(EndT:EndB,EndL) = 1;
        img(EndT:EndB,EndR) = 1;
        
        % Top connecting rectangle segments
        img(EndT,EndL:vL-FlnkSep) = 1;
        img(EndT,vR+FlnkSep:EndR) = 1;
        
        % Bottom connecting rectangle segments
        img(EndB,EndL:vL-FlnkSep) = 1;
        img(EndB,vR+FlnkSep:EndR) = 1;
        
        % Cuboid part -------------------------------------
        DiagLen = 42.9190*ones(1,6);
        DiagAngle = 0.8464; 
        DiagAngle = [DiagAngle*ones(1,3) -DiagAngle*ones(1,3)]+pi;
        
        % Line coordinates
        xCoords = [EndL EndL vL-FlnkSep EndR EndR vR+FlnkSep];	% x-Start            
        yCoords = [EndT EndB EndT EndT EndB EndT];              % y-Start         
        xCoords(2,:) = round(DiagLen.*sin(DiagAngle)+xCoords);	% x-End
        yCoords(2,:) = round(DiagLen.*cos(DiagAngle)+yCoords);	% y-End
        
        % Connect the ends
        img(yCoords(2,1):yCoords(2,2),xCoords(2,1)) = 1;
        img(yCoords(2,4):yCoords(2,5),xCoords(2,4)) = 1;
        img(yCoords(2,1),xCoords(2,1):xCoords(2,3)) = 1;
        img(yCoords(2,4),xCoords(2,6):xCoords(2,4)) = 1;
        
        
        % Line parameters  
        theta = atan2(diff(yCoords),diff(xCoords));
        w1 = sin(theta);
        w2 = -cos(theta);
        w0 = -w1.*xCoords(1,:) -w2.*yCoords(1,:);

        % Parameters of lines perpendicular to lines at the
        % line endpoints
        PerpTheta = theta+pi/2;
        w1Perp = sin(PerpTheta);
        w2Perp = -cos(PerpTheta);


        % Imaging lines
        [X,Y] = meshgrid(1:ImSz(2),1:ImSz(1));
        nLn = length(DiagLen);
        TheDist = inf([ImSz nLn]);
        for iLn = 1:nLn

            % Find the mid-point of the line
            xMid = mean(xCoords(:,iLn));
            yMid = mean(yCoords(:,iLn));

            ValInds = true(ImSz);
            for iEndPt = 1:2

                % Get the final perpendicular line parameter
                w0Perp = -w1Perp(iLn).*(xCoords(iEndPt,iLn)) - w2Perp(iLn).*(yCoords(iEndPt,iLn));

                % Distance of each point from the lines perpendicular to the lines
                PerpDist = w0Perp + w1Perp(iLn)*X + w2Perp(iLn)*Y;

                % Get the sign of the perpendicular distance to the line's mid
                % point
                MidDist = w0Perp + w1Perp(iLn)*xMid + w2Perp(iLn)*yMid;

                ValInds(sign(MidDist)*PerpDist<0) = false;

            end                

            % Distance of each point from the chevron line
            TD = TheDist(:,:,iLn);
            TD(ValInds) = w0(iLn) + w1(iLn)*X(ValInds) + w2(iLn)*Y(ValInds);
            TheDist(:,:,iLn) = TD;

        end

        % Exponential function of chevron line distance
        TheExp = exp(-TheDist.^2/(2*S^2));

        img = max(cat(3,img,TheExp),[],3);
        
        if DeVectorized
            
            plot([EndL EndL],[Mid(1)-MidFlnkHeight Mid(1)+MidFlnkHeight],'w-');
            plot([EndR EndR],[Mid(1)-MidFlnkHeight Mid(1)+MidFlnkHeight],'w-');
            
            plot([EndL vL-FlnkSep],[Mid(1)+MidFlnkHeight Mid(1)+MidFlnkHeight],'w-');
            plot([vR+FlnkSep EndR],[Mid(1)+MidFlnkHeight Mid(1)+MidFlnkHeight],'w-');
            
            plot([EndL vL-FlnkSep],[Mid(1)-MidFlnkHeight Mid(1)-MidFlnkHeight],'w-');
            plot([vR+FlnkSep EndR],[Mid(1)-MidFlnkHeight Mid(1)-MidFlnkHeight],'w-');
            
            plot(xCoords,yCoords,'w-');
            plot([xCoords(2,1) xCoords(2,1)],[yCoords(2,1) yCoords(2,2)],'w-');
            plot([xCoords(2,4) xCoords(2,4)],[yCoords(2,4) yCoords(2,5)],'w-');
            plot([xCoords(2,1) xCoords(2,3)],[yCoords(2,1) yCoords(2,1)],'w-');
            plot([xCoords(2,6) xCoords(2,4)],[yCoords(2,4) yCoords(2,4)],'w-');                                                            
            
        end
        
    case 'scrambled cuboids'                               
        
        % Vertical line coordinates
        vX = [-57 58 -96 97;
             -57 58 -96 97]+Mid(2); vX = vX*clarkeFactor/2;
        vY = [-29 -29 -38 -38;
            0 0 -10 -10]+Mid(1); vY = vY*clarkeFactor/2;
        
        % Horizontal line coordinates
        hX = [-6 7 -12 13 -47 48;
            -63 64 -69 70 -103 104]+Mid(2); hX = hX*clarkeFactor/2;
        hY = [3 3 7 7 -23 -23;
            3 3 7 7 -23 -23]+Mid(1); hY = hY*clarkeFactor/2;
        
        for iV = 1:size(vX,2)
            img(vY(1,iV):vY(2,iV),vX(1,iV)) = 1;
        end
        for iH = 1:size(hX,2)
            img(hY(1,iH),min(hX(:,iH)):max(hX(:,iH))) = 1;
        end
        
        % Diagonal line coordinates
        xCoords = [-6 7 -48 49 -97 98;
            -39 39 -78 78 -63 63]+Mid(2); xCoords = xCoords*clarkeFactor/2;
        yCoords = [-5 -5 -21 -21 -33 -33;
            -33 -33 12 12 -6 -6]+Mid(1);  yCoords = yCoords*clarkeFactor/2;             
        
        
        % Line parameters  
        theta = atan2(diff(yCoords),diff(xCoords));
        w1 = sin(theta);
        w2 = -cos(theta);
        w0 = -w1.*xCoords(1,:) -w2.*yCoords(1,:);

        % Parameters of lines perpendicular to lines at the
        % line endpoints
        PerpTheta = theta+pi/2;
        w1Perp = sin(PerpTheta);
        w2Perp = -cos(PerpTheta);


        % Imaging lines
        [X,Y] = meshgrid(1:ImSz(2),1:ImSz(1));
        nLn = size(xCoords,2);
        TheDist = inf([ImSz nLn]);
        for iLn = 1:nLn

            % Find the mid-point of the line
            xMid = mean(xCoords(:,iLn));
            yMid = mean(yCoords(:,iLn));

            ValInds = true(ImSz);
            for iEndPt = 1:2

                % Get the final perpendicular line parameter
                w0Perp = -w1Perp(iLn).*(xCoords(iEndPt,iLn)) - w2Perp(iLn).*(yCoords(iEndPt,iLn));

                % Distance of each point from the lines perpendicular to the lines
                PerpDist = w0Perp + w1Perp(iLn)*X + w2Perp(iLn)*Y;

                % Get the sign of the perpendicular distance to the line's mid
                % point
                MidDist = w0Perp + w1Perp(iLn)*xMid + w2Perp(iLn)*yMid;

                ValInds(sign(MidDist)*PerpDist<0) = false;

            end                

            % Distance of each point from the chevron line
            TD = TheDist(:,:,iLn);
            TD(ValInds) = w0(iLn) + w1(iLn)*X(ValInds) + w2(iLn)*Y(ValInds);
            TheDist(:,:,iLn) = TD;

        end

        % Exponential function of chevron line distance
        TheExp = exp(-TheDist.^2/(2*S^2));

        img = max(cat(3,img,TheExp),[],3);
        
    case {'|x|','rectangles |x|'}
        
        EndL = vL-FlnkSep-RecWidth;
        EndR = vR+FlnkSep+RecWidth;
        
        % Outer vertical line segments
        img(Mid(1)-MidFlnkHeight:Mid(1)+MidFlnkHeight,EndL) = 1;
        img(Mid(1)-MidFlnkHeight:Mid(1)+MidFlnkHeight,EndR) = 1;                
        
        % X's
        xOffst = round(((vL-FlnkSep)-EndL)/4);
        xX = [EndL+xOffst (vL-FlnkSep)-xOffst];
        yX = [Mid(1)-MidFlnkHeight Mid(1)+MidFlnkHeight];
        Wdth = 0.4;
        s = 1;
        
        % Cartesian coordinate system over the image
        [x,y] = meshgrid(1:ImSz(2),1:ImSz(1));
        
        % Left-side X
        c = (yX(2)-yX(1))/(xX(1)-xX(2));
        w1 = c/(1+c);
        w2 = 1-w1;
        w0 = -yX(1)*w2 - w1*xX(1);               
        
        d = abs(w0 + w1*x + w2*y)/sqrt(w1^2 + w2^2);
        
        img(exp(-d/(2*s^2))>exp(-Wdth*s/(2*s^2))) = 1;
        
        c = (yX(2)-yX(1))/(xX(2)-xX(1));
        w1 = c/(1+c);
        w2 = 1-w1;
        w0 = -yX(2)*w2 - w1*xX(1);               
        
        d = abs(w0 + w1*x + w2*y)/sqrt(w1^2 + w2^2);
        s = 1;
        img(exp(-d/(2*s^2))>exp(-Wdth*s/(2*s^2))) = 1;
        
        % Right-side X
        xX = [(vR+FlnkSep)+xOffst EndR-xOffst];        
        
        c = (yX(2)-yX(1))/(xX(1)-xX(2));
        w1 = c/(1+c);
        w2 = 1-w1;
        w0 = -yX(1)*w2 - w1*xX(1);               
        
        d = abs(w0 + w1*x + w2*y)/sqrt(w1^2 + w2^2);
        s = 1;
        img(exp(-d/(2*s^2))>exp(-Wdth*s/(2*s^2))) = 1;
        
        c = (yX(2)-yX(1))/(xX(2)-xX(1));
        w1 = c/(1+c);
        w2 = 1-w1;
        w0 = -yX(2)*w2 - w1*xX(1);               
        
        d = abs(w0 + w1*x + w2*y)/sqrt(w1^2 + w2^2);
        s = 1;
        img(exp(-d/(2*s^2))>exp(-Wdth*s/(2*s^2))) = 1;
        
        img(y<yX(1)) = 0;
        img(y>yX(2)) = 0;
        
        if strcmpi(StimType,'rectangles |x|')
            
            % Top connecting rectangle segments
            img(Mid(1)+MidFlnkHeight,EndL:vL-FlnkSep) = 1;
            img(Mid(1)+MidFlnkHeight,vR+FlnkSep:EndR) = 1;
            
            % Bottom connecting rectangle segments
            img(Mid(1)-MidFlnkHeight,EndL:vL-FlnkSep) = 1;
            img(Mid(1)-MidFlnkHeight,vR+FlnkSep:EndR) = 1;
        end                           
        
    otherwise
        
        error('Unknown stimulus type.  Available types include: Vernier, Lines, Rectangles, Cuboids and Scrambled Cuboids.');
end










