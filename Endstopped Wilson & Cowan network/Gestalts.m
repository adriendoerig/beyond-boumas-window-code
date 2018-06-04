function [Img,Data] = Gestalts(offset)

% GESTALTS Get the images and data for the Gestalt experiments.
% 
% Usage: [Img,Data] = Gestalts;
% 
% e.g.
% 
% [Img,Data] = Gestalts;
% 
% figure;
% subplot(2,1,1);
% BarH = bar(1:size(Data,2),Data(1,:)); set(BarH,'FaceColor','k');
% hold on;
% errorbar(1:size(Data,2),Data(1,:),Data(2,:),'k.');
% hold off;
% for ind = 1:size(Data,2)
%     subplot(2,size(Data,2),ind+4);
%     imagesc(Img(:,:,ind));
%     axis image;
%     axis off;
% end
% colormap gray;
% 
% Created by Aaron Clarke, EPFL, LPSY, Herzog Lab, September 5, 20014

Data = [900.4875 273.05 1140.775 391.825 279.5875 832.9375; ...
    383.0315883 95.07037087 371.1198759	118.0627349	93.53836195	306.6171577]

Img = zeros(210,450,7);

% Rectangles
Img(:,:,1) = GestaltStim('Vernier',offset);
Img(:,:,2) = GestaltStim('Lines',offset);
Img(:,:,3) = GestaltStim('Rectangles',offset);
Img(:,:,4) = GestaltStim('|x|',offset);
Img(:,:,5) = GestaltStim('Rectangles |x|',offset);
Img(:,:,6) = GestaltStim('cuboids',offset);
Img(:,:,7) = GestaltStim('scrambled cuboids',offset);

