% use this to collect all model outputs and fit a single sigmoid function.
%
% A.Doerig, LPSY, September 2016

clear
close all

stims = {'circles', 'gestalts', 'hexagons', 'irreg1', 'irreg2', ...
    'malaniaEqual', 'malaniaLong', 'malaniaShort', 'octagons' ...
    'pattern2', 'patternIrregular', 'patternStars', 'squares', 'stars'};

makeDataNew;
FitParams = [];
ThreshElev = [];
outputs = [];

%% run the data set used for the psychometric function first
fitting = 1;

for stim = 1:length(stims)
    stims{stim}
    ThreshElev = [ThreshElev, dataStruct.(matlab.lang.makeValidName(stims{stim}))(1,:)];
    outputs = [outputs, aaronModel(stims{stim})];
    stims{stim};
end

FitParams = FitSigmoidArbMaxMin([outputs 1],[ThreshElev 11.5],100*ones(1,length(ThreshElev(1,:))+1)); % i added the vernier alone data here.

figure(1000)
x = min(outputs):0.005:1;
s=SigmoidML(FitParams,x);
plot(x, s)
hold on
scatter([outputs 1], [ThreshElev 11.5]);
hold off
title(['Psychometric function for ALL stimuli'])
saveas(gcf,['results/psychFunction,ALLstimuli.jpg'])

fitting = 0;

%% Now run all conditions with the computed psychometric function
for i = 1:length(stims)
    aaronModel(stims{i});
end


