

mapsShifted_cum = [];
tuningInd_cum = [];
if local
    stimStr = 'DriftingSine'
else
    stimStr = '25degPerSec_25degPerCycle'
end


%cycle through recordings, rois
for ind = 1 : numel( basePaths )
    basePath = basePaths{ind}
    flyInd = flyInds{ind}
    
    for roiInd = roiInds{ind}
        
        stackInds = GetStackInd( basePath, flyInd, stimStr)
        
        functionalChannel = 1;
        plotFlag = 0;
        generateFitData = false;
        upSampleFactor = 1;
        labelStruct = [];
        roiStruct = [];
        
        [imagingStruct, registrationStruct, maskStruct, timingStruct, matchStruct, stimParamsMetadata, stimParams, tuningStruct, roiInfo] = ...
            LoadSession(basePath, flyInd, stackInds, true);
        
        imagingStruct.stackRaw = imagingStruct.stackRaw(:, roiInd, :);
        
        SPARC_plotTCs
        
    end
end


%% calculate mean DSI

%%put array in 3rd dimension, for code that follows
clear mapsShifted_cum_mean_temp
clear mapsShifted_cum_mean
mapsShifted_cum_mean_temp = mean(mapsShifted_cum);
mapsShifted_cum_std_temp = std(mapsShifted_cum);
mapsShifted_cum_std_half = mapsShifted_cum_std_temp ./ 2
for n = 1:size(mapsShifted_cum_mean_temp, 2);
    mapsShifted_cum_mean(1, 1, n) = mapsShifted_cum_mean_temp(n);
end

domain = [0 45 90 135 180 225 270 315]

weightedXCoords = bsxfun( @times, mapsShifted_cum_mean, reshape( cos( domain * pi / 180 ), [1 1 size( domain, 2 )] ) );
weightedYCoords = bsxfun( @times, mapsShifted_cum_mean, reshape( sin( domain * pi / 180 ), [1 1 size( domain, 2 )] ) );
vecSumX = sum( weightedXCoords, 3 );
vecSumY = sum( weightedYCoords, 3 );
[vecSumTheta, vecSumRho] = cart2pol( vecSumX, vecSumY );
vecAvgX = vecSumX ./ sum( mapsShifted_cum_mean, 3 );
vecAvgY = vecSumY ./ sum( mapsShifted_cum_mean, 3 );
[vecAvgTheta, vecAvgRho] = cart2pol( vecAvgX, vecAvgY );

tcParams_new.bestDir = vecAvgTheta * 180 / pi;
tcParams_new.vecAvgRho = vecAvgRho;
tcParams_new.tuningInd = vecSumRho ./ sum( abs( mapsShifted_cum_mean ), 3 ); % circular variance. this is OSI (Bonhoeffer et al. 1995 Euro. J. Neurosci.), modified for signed responses
tcParams_new.sumResp = sum( mapsShifted_cum_mean, 3 ); % the sign of this quantity indicates whether responses are ON or OFF. if this is zero, then it is indeterminate whether the tuning curve is upward or downward facing
tcParams_new.sumAbsResp = sum( abs( mapsShifted_cum_mean ), 3 );
tcParams_new.meanResp = tcParams_new.sumResp / size( mapsShifted_cum_mean, 3 );
tcParams_new.meanAbsResp = tcParams_new.sumAbsResp / size( mapsShifted_cum_mean, 3 );

%%
figure
errorbar(squeeze(mapsShifted_cum_mean), mapsShifted_cum_std_half, colorz, 'LineWidth', 2);
title({[titley ' bouts: ' num2str(bouties)]; ['mean tuning curve' ' DSI: ' num2str(tcParams_new.tuningInd)]})
ylim([0 1.3])
xlim([0 9])

figure
scatter(ones(1, length(tuningInd_cum)), [tuningInd_cum], colorz);
title({[titley ' bouts: ' num2str(bouties)]; 'DSI distribution'})
hold on
scatter(1, (mean(tuningInd_cum)), 'filled', colorz);
ylim([0 1])
xlim([0 2])
hold off

if strcmp(titley, 'FlpOut') || strcmp(titley, 'FlpOut_short')
    mapsShifted_cum_mean_one = mapsShifted_cum_mean
    mapsShifted_cum_std_half_one = mapsShifted_cum_std_half
    tuningInd_cum_one = tuningInd_cum
end

%%
if saveFigures
    SaveOpenFigures ('/scratch/users/wienecke/figures/SPARC_figure_190612/');
end
