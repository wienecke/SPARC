% plotChrimsonAbsoluteValueResponse
% written for analysis of the SPARC chrimson data set
% Yvette Fisher 6/2019

%% Plot each evoked response from this trial and calculate average amplitude
close all;
ephysSettings;

FigHand = figure('Position',[50, 50, 1800, 800]);
set(gcf, 'Color', 'w');
set(gcf,'renderer','Painters')
set(gca,'TickDir','out'); % The only other option is 'in'


% check if scaled voltage exists
if(isfield(data,'scaledVoltage'))
    dataTrace = data.scaledVoltage; % for current clamp traces
else
    dataTrace = data.scaledCurrent;
end 

timeArray = (1  :  length(data.current) ) / settings.sampRate; % seconds

DURATION_TO_PLOT_BEFORE_FLASH = 1;
DURATION_TO_PLOT_PER_EPOCH = 4; % sec

EPOCH_TO_FIND_PEAK = 0.5;
EPOCH_TO_FIND_BASELINE = 1;

STEP_ONSET_MAGNITUDE = 0;
OFFSET_FOR_PLOT = 0;

stimChanges =  diff( stimulus.shutterCommand );
% find index where the shuttle value increases
epochStartInds = find( stimChanges > STEP_ONSET_MAGNITUDE);

traces = [];

preStimAve =[];
postStimMin = [];

for i = 1 : numel(epochStartInds) - 1    
    currEpochStart = epochStartInds(i) - ( DURATION_TO_PLOT_BEFORE_FLASH * settings.sampRate);
    currEpochStimStart = epochStartInds(i);
    currEpochEnd = epochStartInds(i) + ( DURATION_TO_PLOT_PER_EPOCH * settings.sampRate);
    
    currTimeArray = timeArray( currEpochStart : currEpochEnd) ;
    currTimeArray = currTimeArray - currTimeArray(1); % offset to start at t = 0 for plot;
    
    currVoltage = dataTrace( currEpochStart : currEpochEnd);
    currVoltage = currVoltage + (OFFSET_FOR_PLOT * -1*(i - 1) );% offset to see all the traces over eachother
    
    plot(currTimeArray, currVoltage); hold on;
    traces(i, :) = currVoltage;
    box off;
    
    currStim = stimulus.shutterCommand( currEpochStart: currEpochEnd);
    STIM_OFFSET = -95;
    STIM_SCALE = 100;
    plot(currTimeArray, (currStim* STIM_SCALE) + STIM_OFFSET , 'k'); hold on;
    
end

epochStartInd = 1;
epochStimStart = DURATION_TO_PLOT_BEFORE_FLASH * settings.sampRate;
epochEndInd = epochStimStart + EPOCH_TO_FIND_PEAK * settings.sampRate;

meanTrace = mean( traces );
% median filter the trace and plot average over the raw data:
ORDER_MEDFILTER = 400; % parameter than MM uses
filteredMeanTrace = medfilt1(meanTrace, ORDER_MEDFILTER, 'truncate' )'; % Median filtering of the trace

% find peak of the meanTrace after the epoch
preStimMeanTrace = filteredMeanTrace( epochStartInd : epochStimStart);
preStimMeanVoltage = mean( preStimMeanTrace );


postStimMeanTrace = filteredMeanTrace( epochStimStart : epochEndInd);
meanSubtractedPostResp = postStimMeanTrace - preStimMeanVoltage;

[ ~ , absMaxIndex ] = max( abs( meanSubtractedPostResp ) );
meanAmpResponse = meanSubtractedPostResp( absMaxIndex );


% add average IPSP plotp
plot( currTimeArray,filteredMeanTrace , '-k' );

PLOT_SCALE_OFFSET = 10; % mV
 Vm_mean = mean( preStimAve );
ylim([ -55 , -10])
xlabel('s')
ylabel('mV')
title( [ 'Average evoked response amplitude : '   num2str( round( meanAmpResponse ,2) ) 'mV' ] )

meanAmpResponse
preStimMeanVoltage


%%  Strip chart of the peak chrimson response:
tdTomatoPeakAmp = [18.5422000000000;13.3964000000000;13.8769000000000;19.5533000000000];
GFPPeakAmp = [-10.1924000000000;-5.21040000000000;-3.00150000000000];

data{1} = tdTomatoPeakAmp;
data{2} = GFPPeakAmp;

figure;
set(gcf, 'Color', 'w');
set(gca,'TickDir','out'); % The only other option is 'in'

plotSpread( data , 'showMM', 2 );
ylim([ -15,30]);
ylabel('mV')
niceaxes
title('SPARC chrimson response, mean red cross');

