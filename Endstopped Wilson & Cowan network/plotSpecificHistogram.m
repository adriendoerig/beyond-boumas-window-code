% use this to plot a histogram of results for one set of stimuli using another set for the psych function fit.
%
% A.Doerig, LPSY, September 2016

clear
close all

fittingStims = {'squares'};
psychStimName = {'squares'}; % name of the stimuli set you are using for the fitting
histogramStim = {'stars'};

labels = {'vernier','1star','7stars'};
labelsX = 1:length(labels);

makeDataNew;
FitParams = [];
ThreshElev = [];
outputs = [];
histogramStimsOutput = [];
histogramThreshElev = [];

%% run the data set used for the psychometric function first
fitting = 1;

for run = 1:length(fittingStims)
    fittingStims{run}
    outputs = [outputs, aaronModel(fittingStims{run})];
    ThreshElev = [ThreshElev, dataStruct.(matlab.lang.makeValidName(fittingStims{run}))(1,:)];
    fittingStims{run};
end

FitParams = FitSigmoidArbMaxMin([outputs 1],[ThreshElev 11.5],100*ones(1,length(ThreshElev(1,:))+1)); % i added the vernier alone data here.

histogramStimsOutput = SigmoidML(FitParams,1);
histogramStimsOutput = [histogramStimsOutput, SigmoidML(FitParams,aaronModel(histogramStim{1}))];
histogramThreshElev = dataStruct.(matlab.lang.makeValidName(histogramStim{1}))(1,:);
% errorBars = [43 dataStruct.(matlab.lang.makeValidName(histogramStim{1}))(2,:)];
% histogramErrorBars = zeros(1,2*length(errorBars));
% 
% it = 1;
% for i = 1:length(histogramErrorBars)
%     if mod(i,4) == 1
%         histogramErrorBars(i) = errorBars(it);
%         it = it+1;
%     end
% end

figure(1000)

bars = [[110 histogramThreshElev].', histogramStimsOutput.'];  % the first number is the human data vernier threshold
% barsT = bars.';
bar(labelsX,bars);
set(gca,'xticklabel',labels)
ylabel('Threshold')
% hold on; errorbar(1:length(histogramErrorBars),barsT(:).',histogramErrorBars,'k.'); hold off;
legend('humans', 'Wilson & Cowan network')
% title([histogramStim{1}, ' with ', psychStimName{1}, ' psychometric function'])
saveas(gcf,['results/', histogramStim{1}, ' with ', psychStimName{1}, ' psychometric function HISTOGRAM.jpg'])






