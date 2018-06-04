% use this to collect model outputs for a given set of stimuli and plot them
% with a sigmoid function from another set of stimuli.
%
% A.Doerig, LPSY, September 2016

clear
close all

fittingStims = {'malaniaLong'};
psychStimName = {'malaniaLong'}; % name of the stimuli set you are using for the fitting
scatterStim = {'squares'};

makeDataNew;
FitParams = [];
ThreshElev = [];
outputs = [];
scatterStimsOutput = [];
scatterThreshElev = [];

%% run the data set used for the psychometric function first
fitting = 1;

for stim = 1:length(fittingStims)
    fittingStims{stim}
    outputs = [outputs, aaronModel(fittingStims{stim})];
    ThreshElev = [ThreshElev, dataStruct.(matlab.lang.makeValidName(fittingStims{stim}))(1,:)];
    fittingStims{stim};
end

FitParams = FitSigmoidArbMaxMin([outputs 1],[ThreshElev 11.5],100*ones(1,length(ThreshElev(1,:))+1)); % i added the vernier alone data here.

scatterStimsOutput = aaronModel(scatterStim{1});
scatterThreshElev = dataStruct.(matlab.lang.makeValidName(scatterStim{1}))(1,:);

figure(1000)
x = 0:0.005:1;
s=SigmoidML(FitParams,x);
plot(x, s)
hold on
scatter([outputs 1], [ThreshElev 110.6],'b');
scatter([scatterStimsOutput 1], [scatterThreshElev 104.6],'r');
hold off
title([scatterStim{1}, ' with ', psychStimName{1}, ' psychometric function'])
legend([psychStimName{1}, ' psychometric fct'], psychStimName{1}, scatterStim{1})
saveas(gcf,['results/', scatterStim{1}, ' with ', psychStimName{1}, ' psychometric function.jpg'])






