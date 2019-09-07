function [ trialAveragedImagingStruct, trialAveragedTimingStruct ] = ...
    TrialAverage( imagingStruct, registrationStruct, timingStruct, matchStruct, ...
    channel, upSampleFactor, boutIndsToUse, stimulusPeriodsToUse, stimulusPeriod, ...
    shiftSign, showNSeconds, collapseEpochs, customShift, filterPeriod, misOffset )

trialAveragedTimingStruct = timingStruct; % initialize

trialAveragedTimingStruct.functionRunTime = [trialAveragedTimingStruct.functionRunTime; clock];

tic;

%% set up channel

if ~exist( 'channel', 'var' ) || isempty( channel )
    channel = 2;
end

%% set up upSampleFactor

if ~exist( 'upSampleFactor', 'var' ) || isempty( upSampleFactor )
    upSampleFactor = 1;
end

%% set up boutIndsToUse

if ~exist( 'boutIndsToUse', 'var' )
    boutIndsToUse = [];
end

%% set up stimulusPeriodsToUse

if ~exist( 'stimulusPeriodsToUse', 'var' )
    stimulusPeriodsToUse = [];
end

%% set up stimulusPeriodsToUse

if ~exist( 'stimulusPeriod', 'var' )
    stimulusPeriod = [];
end

%% set up shiftSign

if ~exist( 'shiftSign', 'var' )
    shiftSign = '';
end

%% set up showNSeconds

if ~exist( 'showNSeconds', 'var' )
    showNSeconds = []; % by default if showNSeconds is empty then show all periods
end

%% set up collapseEpochs

if ~exist( 'collapseEpochs', 'var' ) || isempty( collapseEpochs )
    collapseEpochs = false;
end

%% set up collapseEpochs

if ~exist( 'customShift', 'var' ) || isempty( customShift )
    customShift = 0;
end


%% set up filterPeriod

if ~exist( 'filterPeriod', 'var' ) || isempty( filterPeriod )
    filterPeriod = 50;
end


%% set up misOffset flag

if ~exist( 'misOffset', 'var' ) || isempty( misOffset )
    misOffset = false;
end


%% realign stack as needed/demanded

if (~imagingStruct.isAligned || ~imagingStruct.isAlignedAbs) && exist( 'registrationStruct', 'var' ) && ~isempty( registrationStruct )
    imagingStruct = AlignImagingStruct(imagingStruct, registrationStruct);
    
end
stack = imagingStruct.stackRaw(:, :, :, channel);

%% preprocess stack to remove slow changes/background


%stack = HighPassFilter_NaN( stack, 50 / (imagingStruct.frameDuration / 1000), 'gaussian', 0 );

stack = HighPassFilter_NaN( stack, filterPeriod / (imagingStruct.frameDuration / 1000), 'gaussian', 0 );


%%

if isempty(imagingStruct.frameDeltaTs)
    windowSize = imagingStruct.frameDuration / upSampleFactor / 1000.0;
else
    windowSize = mean( diff( imagingStruct.frameDeltaTs ) ) / upSampleFactor;
end

[meanStacks, binEdges, elementsPerBin] = MakeFlattenGridStacks_Resample( ...
    stack, ...
    imagingStruct, ...
    timingStruct, ...
    matchStruct, ...
    windowSize, ...
    @(x) nanmean( x, 3 ), ...
    boutIndsToUse, ...
    stimulusPeriodsToUse, ...
    stimulusPeriod, ...
    shiftSign, ...
    showNSeconds, ...
    collapseEpochs, ...
    customShift, ...
    misOffset);
for epochInd = 2 : numel( binEdges )
    binEdges{epochInd} = binEdges{epochInd - 1}(end) + binEdges{epochInd};
end
binEdges(1 : end - 1) = cellfun( @(x) x(1 : end - 1), binEdges(1 : end - 1), 'UniformOutput', false ); % retain the righthand edge of the last epoch

if collapseEpochs
    for nInd = 3:10
        timingStruct.boutEpochInds(timingStruct.boutEpochInds == nInd) = 3; 
    end
    %timingStruct.boutEpochInds((timingStruct.boutEpochInds ~= 0) & (timingStruct.boutEpochInds ~= 1)) = 1;
end

 %if upSampleFactor == 1 && (~exist( 'stimulusPeriodsToUse', 'var' ) || isempty( stimulusPeriodsToUse ))
 %    gridStacks = FormGridStacks( stack, timingStruct, [], -1, [0 0] );
 %    gridMeanStacks = cell( size( gridStacks ) );
 %    for epochInd = 1 : size( meanStacks, 2 ) % size( gridStacks, 2 )
 %        epochStackSize = size ( gridStacks{1, epochInd} );
 %        epochMeanStack = CropArray( meanStacks{epochInd}, epochStackSize );
 %        gridMeanStacks(~cellfun( 'isempty', gridStacks(:, epochInd) ), epochInd) = {epochMeanStack};
 %    end
 %    
 %    gridStacksNoise = cellfun( @( x, y ) x - y, gridStacks, gridMeanStacks, 'UniformOutput', false );
 %    stackNoise = cell2mat( reshape( gridStacksNoise, [1 1 numel(gridStacksNoise)] ) );
 %end

dataFrameEpochInds = cell2mat( cellfun( @(x, y) ones( [1 size( x, 3 )] ) * y, ...
    meanStacks, num2cell( unique( timingStruct.boutEpochInds )' ), 'UniformOutput', false ) );
%dataFrameEpochIndsu = cell2mat( cellfun( @(x, y) ones( [1 size( x, 3 )] ) * y, ...
%    meanyStacks, mavis, 'UniformOutput', false ) );

trialAveragedTimingStruct.dataFrameStartTimes = cell2mat( binEdges )';
trialAveragedTimingStruct.boutStartTimes = cellfun( @(x) x(1), binEdges )';
trialAveragedTimingStruct.boutStartTimes(end + 1) = trialAveragedTimingStruct.dataFrameStartTimes(end);
trialAveragedTimingStruct.boutDurations = diff( trialAveragedTimingStruct.boutStartTimes );
trialAveragedTimingStruct.boutEpochInds = [0 : size( meanStacks, 2 ) - 1]';
trialAveragedTimingStruct.nBouts = size( trialAveragedTimingStruct.boutDurations, 1 );
trialAveragedTimingStruct.nDataFrames = numel( trialAveragedTimingStruct.dataFrameStartTimes ) - 1;

trialAveragedTimingStruct.boutDataFrameInds = cell( size( meanStacks, 2 ), 1 );
trialAveragedTimingStruct.boutInds = cell( size( meanStacks, 2), 1 );
trialAveragedTimingStruct.stimFrameDurations = cell( size( meanStacks, 2 ), 1 );
trialAveragedTimingStruct.stimFrameStartTimes = cell( size( meanStacks, 2 ), 1 );
trialAveragedTimingStruct.stimParam1 = cell( size( meanStacks, 2 ), 1 );
trialAveragedTimingStruct.stimParam2 = cell( size( meanStacks, 2 ), 1 );
trialAveragedTimingStruct.stimParam3 = cell( size( meanStacks, 2 ), 1 );

for epochInd = 1 : size( meanStacks, 2 )
    trialAveragedTimingStruct.boutDataFrameInds{epochInd} = [find( dataFrameEpochInds == epochInd - 1 )];
    
    firstBoutInd = find( timingStruct.boutEpochInds == epochInd - 1, 1, 'first' );
    firstBoutSampleInds = find( timingStruct.boutInds == firstBoutInd );
    
    trialAveragedTimingStruct.boutInds{epochInd} = timingStruct.boutInds(firstBoutSampleInds);

    if ~imagingStruct.isTrialAveraged && ~any( ismember( imagingStruct.functionCalled, 'JoinExperiments_TrialAverage') )
        trialAveragedTimingStruct.stimFrameDurations{epochInd} = timingStruct.stimFrameDurations(firstBoutSampleInds); % stimFrameDurations is normally calculated by taking a diff of stimFrameStartTimes. however, we are massaging the stimFrameStartTimes to align with the dataFrameStartTimes
        trialAveragedTimingStruct.stimFrameStartTimes{epochInd} = timingStruct.stimFrameStartTimes(firstBoutSampleInds);
        trialAveragedTimingStruct.stimFrameStartTimes{epochInd} = ...
            trialAveragedTimingStruct.stimFrameStartTimes{epochInd} - trialAveragedTimingStruct.stimFrameStartTimes{epochInd}(1) + ...
            trialAveragedTimingStruct.boutStartTimes(epochInd);
        trialAveragedTimingStruct.stimParam1{epochInd} = timingStruct.stimParam1(firstBoutSampleInds);
        trialAveragedTimingStruct.stimParam2{epochInd} = timingStruct.stimParam2(firstBoutSampleInds);
        trialAveragedTimingStruct.stimParam3{epochInd} = timingStruct.stimParam3(firstBoutSampleInds);
    end
end

trialAveragedTimingStruct.stimFrameDurations = cell2mat( trialAveragedTimingStruct.stimFrameDurations );
trialAveragedTimingStruct.stimFrameStartTimes = cell2mat( trialAveragedTimingStruct.stimFrameStartTimes );
trialAveragedTimingStruct.stimParam1 = cell2mat( trialAveragedTimingStruct.stimParam1 );
trialAveragedTimingStruct.stimParam2 = cell2mat( trialAveragedTimingStruct.stimParam2 );
trialAveragedTimingStruct.stimParam3 = cell2mat( trialAveragedTimingStruct.stimParam3 );

trialAveragedTimingStruct.boutDataFrameIndsCropped = trialAveragedTimingStruct.boutDataFrameInds;
[~, trialAveragedTimingStruct.boutStartSampleInds, trialAveragedTimingStruct.boutInds] = ...
    unique( cell2mat( trialAveragedTimingStruct.boutInds ), 'stable' );
trialAveragedTimingStruct.boutStartSampleInds(end + 1) = numel( trialAveragedTimingStruct.stimFrameDurations ) + 1;
trialAveragedTimingStruct.gridDataFrameInds = trialAveragedTimingStruct.boutDataFrameInds';
trialAveragedTimingStruct.gridDataFrameIndsCropped = trialAveragedTimingStruct.gridDataFrameInds;
trialAveragedTimingStruct.isTrialAveraged = true;
trialAveragedTimingStruct.nStimFrames = size( trialAveragedTimingStruct.stimFrameDurations, 1 );
trialAveragedTimingStruct.stimEpochInds = trialAveragedTimingStruct.boutInds - 1;
trialAveragedTimingStruct.stimFrameInds = [1 : trialAveragedTimingStruct.nStimFrames]';

trialAveragedTimingStruct.boutIndsUnique = [];
trialAveragedTimingStruct.dataFrameDurations = [];
trialAveragedTimingStruct.dataFrameInds = [];
trialAveragedTimingStruct.dataFrameIndsUnique = [];
trialAveragedTimingStruct.dataFrameStartSampleInds = [];
trialAveragedTimingStruct.orphanDataFrameInds = [];
trialAveragedTimingStruct.stimulusOutputRaw = '';
trialAveragedTimingStruct.boutIndsToUse = boutIndsToUse;
trialAveragedTimingStruct.stimulusPeriod = stimulusPeriod;
trialAveragedTimingStruct.stimulusPeriodsToUse = stimulusPeriodsToUse;

functionElapsedTime = toc;
[ST, ~] = dbstack;
disp( [ ST(1).name ': ' num2str( functionElapsedTime ) ' s'] );
try trialAveragedTimingStruct.functionCalled{end+1} = ST(1).name;
    catch trialAveragedTimingStruct.functionCalled = {trialAveragedTimingStruct.functionCalled ST(1).name}; 
end
trialAveragedTimingStruct.functionElapsedTime = [trialAveragedTimingStruct.functionElapsedTime ; functionElapsedTime];
trialAveragedTimingStruct = orderfields( trialAveragedTimingStruct );

trialAveragedImagingStruct = imagingStruct;
trialAveragedImagingStruct.stackRaw = cell2mat( reshape( meanStacks, [1 1 numel( meanStacks )] ) );
 %if upSampleFactor == 1 && (~exist( 'stimulusPeriodsToUse', 'var' ) || isempty( stimulusPeriodsToUse ))
 %    trialAveragedImagingStruct.gridStacksNoise = gridStacksNoise;
 %    trialAveragedImagingStruct.stackNoise = stackNoise;
 %end
trialAveragedImagingStruct.elementsPerBin = elementsPerBin;
trialAveragedImagingStruct.minEltsPerBin = min( cellfun( @min, elementsPerBin ) );
trialAveragedImagingStruct.maxEltsPerBin = max( cellfun( @max, elementsPerBin ) );
trialAveragedImagingStruct.minEltsPerBin_nonBoundary = min( cellfun( @(x) min( x(2 : end - 1) ), elementsPerBin ) );
trialAveragedImagingStruct.maxEltsPerBin_nonBoundary = max( cellfun( @(x) max( x(2 : end - 1) ), elementsPerBin ) );
trialAveragedImagingStruct.minEltsPerBin_nonBlank = min( cellfun( @min, elementsPerBin(2 : end) ) );
trialAveragedImagingStruct.maxEltsPerBin_nonBlank = max( cellfun( @max, elementsPerBin(2 : end) ) );
trialAveragedImagingStruct.frameDeltaTs = [];
trialAveragedImagingStruct.metadataRaw = '';
trialAveragedImagingStruct.nFrames = size( trialAveragedImagingStruct.stackRaw, 3 );
trialAveragedImagingStruct.physicalDims(3) = trialAveragedImagingStruct.nFrames * trialAveragedImagingStruct.frameDuration / 1000;
trialAveragedImagingStruct.pixelDims(3) = trialAveragedImagingStruct.nFrames;
% trialAveragedImagingStruct.stackAquisitionTime = [];
% trialAveragedImagingStruct.stackInd = [];
% trialAveragedImagingStruct.stackName = '';
trialAveragedImagingStruct.isTrialAveraged = true;
trialAveragedImagingStruct.isAligned = true;
trialAveragedImagingStruct.isAlignedAbs = true;
if size( imagingStruct.stackRaw, 4 ) ~= 1
    trialAveragedImagingStruct.channel = channel;
else
    if isfield( imagingStruct, 'channel' ) && ~isempty( imagingStruct.channel )
        trialAveragedImagingStruct.channel = imagingStruct.channel;
    else
        trialAveragedImagingStruct.channel = [];
    end
end

try trialAveragedImagingStruct.functionCalled{end+1} = ST(1).name;
    catch trialAveragedImagingStruct.functionCalled = {trialAveragedImagingStruct.functionCalled ST(1).name}; 
end
trialAveragedImagingStruct = orderfields( trialAveragedImagingStruct );

end
