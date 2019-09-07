function [roiTCsNormalized, gridStacks, epochStacks] = PlotRoiTimeCourses( ...
    imagingStruct, timingStruct, registrationStruct, roiStruct, ...
    channel, normalizationMethod, roiIndsHighlight, roiLabels, ...
    plotCustomStack, customStack, titley, stimStr, ...
    plotDiff, barDiff, GECODir, roiIndly)

if ~exist( 'plotDiff', 'var' ) || isempty( plotDiff )
    plotDiff = false;
end
if ~exist( 'plotCustomStack', 'var' ) || isempty( plotCustomStack )
    plotCustomStack = false;
end
if ~exist( 'titley', 'var' ) || isempty( titley )
    titley = '';
end
if ~exist( 'stimStr', 'var' ) || isempty( stimStr )
    stimStr = '';
end

if ~exist( 'barDiff', 'var' ) || isempty( barDiff )
    barDiff = {};
end

if ~exist( 'GECODir', 'var' ) || isempty( GECODir )
    GECODir = [];
end

if ~exist( 'roiIndly', 'var' ) || isempty( roiIndly )
    roiIndly = [];
end

if isfield( imagingStruct, 'roiStruct' )
    roiStruct = imagingStruct.roiStruct;
    roiTCs = imagingStruct.stackRaw(:, :, :, channel);
else
    if (~imagingStruct.isAligned || ~imagingStruct.isAlignedAbs) && exist( 'registrationStruct', 'var' ) && ~isempty( registrationStruct )
        imagingStruct = AlignImagingStruct(imagingStruct, registrationStruct);
    end
    stack = imagingStruct.stackRaw(:, :, :, channel);
    
    roiTCs = ApplyROIsToStack_TimeCourse( stack, roiStruct);
    roiTCs = reshape( roiTCs, [1 size( roiTCs )] );
end

if ~imagingStruct.isTrialAveraged
    roiTCs = HighPassFilter_NaN( roiTCs, 50 / (imagingStruct.frameDuration / 1000), 'gaussian', 0 );
    %     roiTCs = HighPassFilter_NaN( roiTCs, 50 / (imagingStruct.frameDuration / 1000), 'gaussian', 0 ); % kludge for filtering to correct fast bleaching. included here but not in other functions because this function is primarily for plotting
end

if isfield( imagingStruct, 'stackRawUnits' ) && ...
        ~isempty( imagingStruct.stackRawUnits ) && ...
        ~strcmp( imagingStruct.stackRawUnits, 'F' )
    if exist( 'normalizationMethod', 'var' )
        disp( ['imagingStruct.stackRawUnits already ' imagingStruct.stackRawUnits '.'] );
        disp( ['No further normalization performed despite normalizationMethod specified as ' normalizationMethod '.'] );
    end
    normalizationMethod = imagingStruct.stackRawUnits;
    roiTCsNormalized = roiTCs;
else
    if ~exist( 'normalizationMethod', 'var' ) || isempty( normalizationMethod )
        normalizationMethod = 'dF_F0';
    end
    roiTCsNormalized = NormalizeStack( roiTCs, timingStruct, normalizationMethod );
    nPredictions = 2;
    if strcmp( timingStruct.functionCalled(end), 'JoinExperiments_TrialAveragey' ) && ...
            max( unique( [timingStruct.boutDataFrameInds{:}] ) ) ~= numel( unique( [timingStruct.boutDataFrameInds{:}] ) ) % check that there is a splice and splice is nontrivial
        if numel( unique( timingStruct.boutEpochInds)) == nPredictions + 1 % if only predictions + 1 control
            [num2str( nPredictions ) ' predictions + 1 blank control']
            roiTCsNormalized(:, :, [timingStruct.boutDataFrameInds{timingStruct.boutEpochInds >= 1}]) = ...
                2 .*  roiTCsNormalized(:, :, [timingStruct.boutDataFrameInds{timingStruct.boutEpochInds >= 1}]);
        elseif numel( unique( timingStruct.boutEpochInds)) == 2 * nPredictions + 1 % otherwise check if there are equal amounts of data and predictions
            roiTCsNormalized(:, :, [timingStruct.boutDataFrameInds{timingStruct.boutEpochInds >= nPredictions + 1}]) = ...
                2 .*  roiTCsNormalized(:, :, [timingStruct.boutDataFrameInds{timingStruct.boutEpochInds >= nPredictions + 1}]);
            
        end
    end
end

if ~exist( 'roiIndsHighlight', 'var' )
    roiIndsHighlight = [];
end

roiLabels = RoundFloat( roiLabels, 2 );
% sensor = 'ASAP2f'; % for hard coding in the event of old or nonconforming
% file names
sensor = GetSensor( imagingStruct, channel );

if HasBlankEpoch( timingStruct ) && ~imagingStruct.isTrialAveraged
    spacerFrames = size( timingStruct.gridDataFrameIndsCropped{1,1}, 2 ); % number of interleaved blank frames
else
    spacerFrames = 0; % no blank
end

[~, gridStacks] = MakeGridStacks( roiTCsNormalized, timingStruct, [], -1 * spacerFrames );
epochStacks = FlattenGridStacks( gridStacks );

if plotDiff == true
    epochStacks = barDiff;
    gridStacks = barDiff;
end
if plotCustomStack == true
    epochStacks = customStack;
    gridStacks = customStack;
    roiNumbe = size(customStack{1}, 2)
    roiTCs(:,roiNumbe,:) = 0;
    roiTCsNormalized(:,roiNumbe,:) = 0;
end


domain = StimulusDomain( timingStruct );

panelRows = 6;
panelCols = 4;

% set up file paths and names
patty = path;
if ~isempty(regexp( patty, '/Users/wienecke/Documents/GitHub/' ));
    imagingStruct(1).experimentParentDir = '~/Documents/ambrose/leprechaun/'
end
if ~isempty(regexp( patty, '/home/users/wienecke' ));
    imagingStruct(1).experimentParentDir = '/scratch/users/wienecke/leprechaun/'
end

roiFileDir = [regexprep( imagingStruct(1).experimentParentDir, ['([^\' filesep ']+)\' filesep '$'], ['$1ROI\' filesep] ) ...
    imagingStruct(1).experimentDir filesep];
if ~isdir( roiFileDir ) mkdir( roiFileDir ); end

% save the roiStruct just in case
roiStructFileNameRoot = ['roiStruct_' roiStruct.uniqueID];
% if exist( [roiFileDir roiStructFileNameRoot '.mat'], 'file' ) ~= 2
%     save( [roiFileDir roiStructFileNameRoot '.mat'], 'roiStruct', '-v7.3', '-mat' );
% end

plotFileNameRoot = [imagingStruct(1).fileNameRoot imagingStruct(1).fileNameSuffix];
plotFilePathRoot = [roiFileDir roiStructFileNameRoot '_' plotFileNameRoot '_' num2str( imagingStruct.stackInd )];

for pageInd = 1 : ceil( size( roiTCs, 2 ) / 4 )
    
    pagedRoiInds = [1 + 4 * (pageInd - 1) : min( size( roiTCs, 2 ), 4 + 4 * (pageInd - 1) )];
    
    figureHandle = figure;
    set( gcf, 'Renderer', 'painters' );
    set( gcf, 'PaperPositionMode', 'manual' );
    set( gcf, 'PaperUnits', 'inches' );
    set( gcf, 'PaperSize', [8.5 11] );
    set( gcf, 'PaperPosition', [0 0.5 8.5 11] );
    
    %%
    
    subplot( panelRows, panelCols, [3 : 4] );
    hold on;
    
    if plotCustomStack == false
        %text( 0, 9, ' ', 'FontSize', 8, 'Interpreter', 'none' );
        text( 0, 8, imagingStruct.fileName, 'FontSize', 8, 'Interpreter', 'none' );
        text( 0, 7, timingStruct.stimParamFileName, 'FontSize', 8, 'Interpreter', 'none' );
        text( 0, 6, ['stack: ' num2str(imagingStruct.stackInd)], 'FontSize', 8, 'Interpreter', 'none' );
        text( 0, 5, ['sensor: ' sensor], 'FontSize', 8, 'Interpreter', 'none' );
        if imagingStruct.isTrialAveraged
            text( 0, 4, ['bouts: ' num2str( timingStruct.boutIndsToUse )], 'FontSize', 8, 'Interpreter', 'none' );
            text( 0, 3, ['period: ' num2str( timingStruct.stimulusPeriod ) ' sec'], 'FontSize', 8, 'Interpreter', 'none' );
            text( 0, 2, ['periods: ' num2str( timingStruct.stimulusPeriodsToUse )], 'FontSize', 8, 'Interpreter', 'none' );
        end
        text( 0, 1, ['roiStruct uniqueID: ' roiStruct.uniqueID], 'FontSize', 8, 'Interpreter', 'none' );
        text( 0, 9, ['GECO bestDir: ' num2str(GECODir) 'roi' num2str(roiIndly)], 'FontSize', 8, 'Interpreter', 'none' );
    end
    
    if plotCustomStack == true
        text( 0, 8, stimStr, 'FontSize', 8, 'Interpreter', 'none' );
        text( 0, 9, titley, 'FontSize', 8, 'Interpreter', 'none' );
    end
    
    ylim( [0 6] );
    axis off;
    hold off;
    
    %% ROI overlay
    
    if plotCustomStack == false
        tmpRoiStruct = roiStruct;
        [tmpRoiStruct.bwMask, tmpRoiStruct.bwLabel, tmpRoiStruct.nROIs, tmpRoiStruct.bwMaskStack] = ...
            UpdateROIData_Subsample( roiStruct.bwMaskStack, pagedRoiInds );
        
        subplot( panelRows, panelCols, [1 : 2] );
        if isfield( imagingStruct, 'roiStruct' )
            labelData = roiStruct.bwLabel / max( roiStruct.bwLabel(:) );
            labelData(labelData == 0) = NaN;
            PlotROIOverlay( labelData, tmpRoiStruct, 'hsv', true );
        else
            PlotROIOverlay( imagingStruct, tmpRoiStruct, 'rgb', true );
        end
    end
    
    %% unbinned timecourse
    
    plotScaleFactor = max( range( roiTCsNormalized(1, pagedRoiInds, :), 3 ) )
    if ~isfinite( plotScaleFactor ) || plotScaleFactor == 0
        plotScaleFactor = 1;
    end
    
    panelInds = [5 : 8];
    subplot( panelRows, panelCols, panelInds );
    hold on;
    
    for boutInd = 1 : timingStruct.nBouts
        if (HasBlankEpoch( timingStruct ) && timingStruct.boutEpochInds(boutInd) ~= 0)
            fill( [timingStruct.boutDataFrameInds{boutInd}(1) timingStruct.boutDataFrameInds{boutInd}(1) timingStruct.boutDataFrameInds{boutInd}(end) timingStruct.boutDataFrameInds{boutInd}(end)], ...
                [1.5 -4 -4 1.5] * plotScaleFactor, ...
                [0.75 0.75 0.75], ...
                'EdgeColor', 'none');
            %             keyboard
            text( timingStruct.boutDataFrameInds{boutInd}(1), plotScaleFactor, ...
                num2str( domain(:, timingStruct.boutEpochInds(boutInd)) ), 'Color', [0 0 0], 'FontSize', 8 );
        elseif (~HasBlankEpoch( timingStruct ) && mod( boutInd, 2 ) == 1)
            fill( [timingStruct.boutDataFrameInds{boutInd}(1) timingStruct.boutDataFrameInds{boutInd}(1) timingStruct.boutDataFrameInds{boutInd}(end) timingStruct.boutDataFrameInds{boutInd}(end)], ...
                [1.5 -4 -4 1.5] * plotScaleFactor, ...
                [0.75 0.75 0.75], ...
                'EdgeColor', 'none');
            text( timingStruct.boutDataFrameInds{boutInd}(1), plotScaleFactor, ...
                num2str( domain(:, timingStruct.boutEpochInds(boutInd) + 1) ), 'Color', [0 0 0], 'FontSize', 8 );
        end
    end
    
    % shade orphan frames red
    if ~isempty( timingStruct.orphanDataFrameInds )
        fill( [timingStruct.orphanDataFrameInds(1) timingStruct.orphanDataFrameInds(1) timingStruct.orphanDataFrameInds(end) timingStruct.orphanDataFrameInds(end)], ...
            [1.5 -4 -4 1.5] * plotScaleFactor, ...
            [0.75 0 0], ...
            'EdgeColor', 'none');
    end
    
    % plot the data as dF/F0. the first trace has dF/F0 = 0 at a y value of 1,
    % the second trace has dF/F0 = 0 at a y value of 0
    for roiInd = pagedRoiInds
        plot( squeeze( roiTCsNormalized(1, roiInd, :) ) - plotScaleFactor * (roiInd - 4 * (pageInd - 1) - 1), 'k-', 'LineWidth', 0.5 );
        
        if ~isempty( strfind( roiIndsHighlight, roiInd ) )
            lineColor = [0.75 0 0.75];
        else
            lineColor = [0.75 0.75 0.75];
        end
        line( [1 size( roiTCsNormalized, 3 )], ones( 1, 2 ) * -1 * plotScaleFactor * (roiInd - 4 * (pageInd - 1) - 1), 'LineStyle', '-', 'Color', lineColor, 'LineWidth', 0.5 );
        text( spacerFrames / 2, -plotScaleFactor * (roiInd - 4 * (pageInd - 1) - 1) + plotScaleFactor / 2, ...
            num2str( roiInd ), 'Color', [0 0 0], 'FontSize', 8 );
    end
    
    fill( [1 1 size(roiTCs, 3) / 1e3 size(roiTCs, 3) / 1e3] + 3 * size(roiTCs, 3) / 1e3, [0 -1 -1 0] + 1.25 * plotScaleFactor, [0 0 0] ); % scale bar
    
    xlim( [1 size(roiTCs, 3)] );
    ylim( [-3.5 1.5] * plotScaleFactor );
    set( gca, 'YTick', [] );
    xlabel( 'data frame index' );
    ylabel( normalizationMethod, 'Interpreter', 'none' );
    hold off;
    
    %% plot condition-sorted single trials and trial averages
    
    if size( domain, 1 ) == 1 || ...
            (size( domain, 1 ) == 2 && ...
            (numel( unique( domain(1, :) ) ) == numel( domain(1, :) ) && ...
            numel( unique( domain(2, :) ) ) == numel( domain(2, :) ))) % if the domain is one-dimensional
        % trial-averaged timecourse
        if imagingStruct.isTrialAveraged
            panelInds = [5 : 20];
        else
            panelInds = [9 : 20];
        end
        subplot( panelRows, panelCols, panelInds );
        hold on;
        
        % shade the stimulus bouts that are not blank
        nonBlankEpochInds = unique( timingStruct.boutEpochInds )';
        if HasBlankEpoch( timingStruct )
            nonBlankEpochInds = nonBlankEpochInds(2 : end);
        end
        frameStart = 1;
        for epochInd = nonBlankEpochInds
            nFrames = size( timingStruct.gridDataFrameIndsCropped{1, epochInd + 1}, 2 );
            if numel( unique( timingStruct.stimParamPerSampleData{epochInd + 1, 1} ) ) == 2 % if it's an ON-OFF stimulus, e.g. LocalCircle
                onFrames = max( 2, ...
                    numel( find( timingStruct.stimParamPerSampleData{epochInd + 1, 1} == timingStruct.stimParamPerSampleData{epochInd + 1, 1}(1) ) ) / ...
                    numel( timingStruct.stimParamPerSampleData{epochInd + 1, 1} ) * ...
                    nFrames );
            else
                onFrames = nFrames;
            end
            
            if ~imagingStruct.isTrialAveraged
                fill( [frameStart frameStart frameStart + onFrames - 1 frameStart + onFrames - 1], ...
                    [1.5 -4 -4 1.5] * plotScaleFactor, ...
                    [0.75 0.75 0.75], ...
                    'EdgeColor', 'none'); % stimulus followed by blank with stimulus in the shaded area and blank the white area
            else
                line( [frameStart frameStart], [-4 1.5] * plotScaleFactor, 'LineStyle', '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5 ); % start line is solid
                line( [frameStart + onFrames - 1 frameStart + onFrames - 1], [-4 1.5] * plotScaleFactor, 'LineStyle', '--', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5 ); % end/half-way line is dashed
            end
            textColor = [0 0 0];
            if HasBlankEpoch( timingStruct )
                text( frameStart, plotScaleFactor, ...
                    num2str( domain(:, epochInd) ), 'Color', textColor, 'FontSize', 8 );
            else
                text( frameStart, plotScaleFactor, ...
                    num2str( domain(:, epochInd + 1) ), 'Color', textColor, 'FontSize', 8 );
            end
            
            frameStart = frameStart + nFrames + spacerFrames;
        end
        
        for roiInd = pagedRoiInds
            frameStart = 1;
            for epochInd = 1 : size( epochStacks, 2 )
                nFrames = size( epochStacks{epochInd}, 3 );
                epochSingleTrials = cat( 1, gridStacks{:, epochInd} );
                plot( [frameStart : frameStart + nFrames - 1], squeeze( epochSingleTrials(:, roiInd, :) ) - plotScaleFactor * (roiInd - 4 * (pageInd - 1) - 1), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1 );
                h(roiInd, epochInd) = ...
                    plot( [frameStart : frameStart + nFrames - 1], squeeze( epochStacks{epochInd}(1, roiInd, :) ) - plotScaleFactor * (roiInd - 4 * (pageInd - 1) - 1), '-', 'Color', [0 0 0], 'LineWidth', 1 );
                
                if ~isempty( strfind( roiIndsHighlight, roiInd ) )
                    lineColor = [0.75 0 0.75];
                else
                    lineColor = [0.75 0.75 0.75];
                end
                line( [frameStart frameStart + nFrames - 1], ones( 1, 2 ) * -1 * plotScaleFactor * (roiInd - 4 * (pageInd - 1) - 1), 'LineStyle', '-', 'Color', lineColor, 'LineWidth', 0.5 );
                
                frameStart = frameStart + nFrames;
            end
            
            text( 0, -plotScaleFactor * (roiInd - 4 * (pageInd - 1) - 1) + plotScaleFactor / 2, ...
                num2str( roiInd ), 'Color', [0 0 0], 'FontSize', 8 );
        end
        
        for roiInd = pagedRoiInds
            for epochInd = 1 : size( epochStacks, 2 )
                uistack( h(roiInd, epochInd), 'top' );
            end
        end
        
        fill( [1 1 frameStart / 1e3 frameStart / 1e3] + 3 * frameStart / 1e3, [0 -1 -1 0] + 1.25 * plotScaleFactor, [0 0 0] ); % scale bar
        
        xlim( [1 frameStart] );
        ylim( [-3.5 1.5] * plotScaleFactor );
        set( gca, 'YTick', [] );
        xlabel( 'data frame index' );
        ylabel( normalizationMethod, 'Interpreter', 'none' );
        hold off;
        
        XTPlotInd = 1;
        for roiInd = pagedRoiInds
            ax1 = subplot( panelRows, panelCols, 20 + XTPlotInd );
            
            maxFramesPerBout = max( cell2mat( ...
                cellfun( @size, epochStacks, repmat( {3}, size( epochStacks ) ), ...
                'UniformOutput', false ) ) );
            XTReceptiveField = zeros( size( epochStacks, 2 ), maxFramesPerBout );
            for epochInd = 1 : size( epochStacks, 2 )
                XTReceptiveField(epochInd, 1 : numel( epochStacks{1, epochInd}(1, roiInd, :) )) = squeeze( epochStacks{1, epochInd}(1, roiInd, :) );
            end
            
            imagesc( XTReceptiveField, ZeroCenteredBounds( XTReceptiveField, [0 100] ) );
            
            colormap( ColorMap_RedWhiteBlue( 129 ) );
            %             colormap( cool );
            %             colorbar;
            axis square;
            if exist( 'roiLabels', 'var' ) && ~isempty( roiLabels )
                if size(roiLabels,2) == 1
                    textString = ['ROI ' num2str( roiInd ) ' (' num2str( roiLabels(roiInd, :), '%0.2f %0.2f' ) ')'];
                    xlabel( textString, 'FontSize', 8 );
                end
                if size(roiLabels,2) == 2
                    textString = ['ROI ' num2str( roiInd ) ' (' num2str( roiLabels(roiInd, :), '%0.2f %0.2f' ) ')'];
                    xlabel( textString, 'FontSize', 8 );
                end
                if size(roiLabels,2) > 2
                    textString = ['ROI ' num2str( roiInd ) ' (' num2str( roiLabels(roiInd, 1:2), '%0.2f %0.2f' ) ')'];
                    textString2 = [' (' num2str( roiLabels(roiInd, 3:4), '%0.2f %0.2f' ) ')'];
                    textString3 = [' (' num2str( roiLabels(roiInd, 5), '%0.2f %0.2f' ) ')'];
                    xlabel( {textString; textString2; textString3}, 'FontSize', 8 );
                end
            else
                textString = ['ROI ' num2str( roiInd )];
                xlabel( textString, 'FontSize', 8 );
            end
            %             set( ax1, 'XAxisLocation', 'bottom', ...
            %                 'YAxisLocation', 'left', ...
            %                 'Color', 'none', ...
            %                 'XColor', 'k', 'YColor', 'k', ...
            %                 'YTickLabel', num2str( domain( get( gca, 'YTick' ) )' ), ...
            %                 'FontSize', 8 );
            set( ax1, 'XAxisLocation', 'bottom', ...
                'YAxisLocation', 'left', ...
                'Color', 'none', ...
                'XColor', 'k', 'YColor', 'k', ...
                'XTickLabel', '', ...
                'YTickLabel', '' );
            
            %             ax2 = axes( 'Position', get( ax1, 'Position' ), ...
            %                 'XAxisLocation', 'bottom', ...
            %                 'YAxisLocation', 'left', ...
            %                 'Color', 'none', ...
            %                 'XColor', 'k', 'YColor', 'k', ...
            %                 'FontSize', 8 );
            %             hold on;
            %
            %             for epochInd = 1 : size( dF_F0Mean, 2 )
            %                 plot( dF_F0Mean{roiInd, epochInd, :}, ':', 'Color', [0 0 0], 'LineWidth', 0.5, 'Parent', ax2 )
            %             end
            %
            %             xlim( [1 maxFramesPerBout] );
            %             ylim( [min( min( XTReceptiveField ) ) - 0.05 max( max( XTReceptiveField ) ) + 0.05] );
            %             xlabel( 'data frame index' );
            %             ylabel( [normalizationMethod ', ROI ' num2str( roiInd )], 'Interpreter', 'none' );
            %             hold off;
            
            XTPlotInd = XTPlotInd + 1;
        end
    elseif size( domain, 1 ) == 2 && ...
            ~(numel( unique( domain(1, :) ) ) == numel( domain(1, :) ) && ...
            numel( unique( domain(2, :) ) ) == numel( domain(2, :) ))
        plotLocations = { ...
            [9 10 13 14], ...
            [11 12 15 16], ...
            [17 18 21 22], ...
            [19 20 23 24] };
        [uniqueDomainX, ~, uniqueDomainXInd] = unique( domain(1, :) );
        [~, ~, uniqueDomainYInd] = unique( domain(2, :) );
        maxFramesPerBout = max( max( cell2mat( ...
            cellfun( @size, timingStruct.gridDataFrameIndsCropped, repmat( {2}, [size( timingStruct.gridDataFrameIndsCropped )] ), ...
            'UniformOutput', false ) ) ) );
        
        XYPlotInd = 1;
        for roiInd = pagedRoiInds
            subplot( panelRows, panelCols, plotLocations{XYPlotInd} );
            hold on
            
            % shade the stimulus bouts that are not blank
            nonBlankEpochInds = unique( timingStruct.boutEpochInds )';
            if HasBlankEpoch( timingStruct )
                nonBlankEpochInds = nonBlankEpochInds(2 : end);
            end
            for epochInd = nonBlankEpochInds
                if HasBlankEpoch( timingStruct )
                    xLoc = maxFramesPerBout * (uniqueDomainXInd(epochInd) - 1);
                    yLoc = plotScaleFactor * (uniqueDomainYInd(epochInd) - 1) + 1;
                else
                    xLoc = maxFramesPerBout * (uniqueDomainXInd(epochInd + 1) - 1);
                    yLoc = plotScaleFactor * (uniqueDomainYInd(epochInd + 1) - 1) + 1;
                end
                nFrames = size( timingStruct.gridDataFrameIndsCropped{1, epochInd + 1}, 2 );
                if mod( xLoc / maxFramesPerBout + 1, 2 ) == 1
                    fill( [xLoc xLoc xLoc + nFrames - 1 xLoc + nFrames - 1], ...
                        [yLoc - 0.5 * plotScaleFactor yLoc + 0.5 * plotScaleFactor yLoc + 0.5 * plotScaleFactor yLoc - 0.5 * plotScaleFactor], ...
                        [0.75 0.75 0.75], ...
                        'EdgeColor', 'none'); % stimulus followed by blank with stimulus in the shaded area and blank the white area
                    if HasBlankEpoch( timingStruct )
                        text( xLoc, yLoc, ...
                            num2str( domain(:, epochInd) ), 'Color', [1 1 1], 'FontSize', 8 );
                    else
                        text( xLoc, yLoc, ...
                            num2str( domain(:, epochInd + 1) ), 'Color', [1 1 1], 'FontSize', 8 );
                    end
                else
                    line( [xLoc xLoc], [yLoc - 0.5 * plotScaleFactor yLoc + 0.5 * plotScaleFactor], 'LineStyle', ':', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5 );
                    line( [xLoc + nFrames xLoc + nFrames], [yLoc - 0.5 * plotScaleFactor yLoc + 0.5 * plotScaleFactor], 'LineStyle', ':', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5 );
                    if HasBlankEpoch( timingStruct )
                        text( xLoc, yLoc, ...
                            num2str( domain(:, epochInd) ), 'Color', [0 0 0], 'FontSize', 8 );
                    else
                        text( xLoc, yLoc, ...
                            num2str( domain(:, epochInd + 1) ), 'Color', [0 0 0], 'FontSize', 8 );
                    end
                end
            end
            
            if ~isempty( strfind( roiIndsHighlight, roiInd ) )
                lineColor = [1 0 1];
            else
                lineColor = [0 0 0];
            end
            for epochInd = 1 : size( epochStacks, 2 )
                xLoc = maxFramesPerBout * (uniqueDomainXInd(epochInd) - 1);
                yLoc = plotScaleFactor * (uniqueDomainYInd(epochInd) - 1) + 1;
                nFrames = size( epochStacks{epochInd}, 3 );
                epochSingleTrials = cat( 1, gridStacks{:, epochInd} );
                plot( [xLoc : xLoc + nFrames - 1], squeeze( epochSingleTrials(:, roiInd, :) ) + yLoc, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5 );
                plot( [xLoc : xLoc + nFrames - 1], squeeze( epochStacks{epochInd}(1, roiInd, :) ) + yLoc, '-', 'Color', lineColor, 'LineWidth', 0.5 );
            end
            
            text( 5, plotScaleFactor * numel( uniqueDomainX ) + 1, ...
                ['ROI ' num2str( roiInd )], 'Color', [0 0 0], 'FontSize', 8 );
            
            xlim( [0 maxFramesPerBout * numel( uniqueDomainX )] );
            ylim( [0.5 plotScaleFactor * numel( uniqueDomainX ) + 1] );
            hold off;
            
            XYPlotInd = XYPlotInd + 1;
        end
    end
    
    %     saveas( figureHandle, [plotFilePathRoot '_page' num2str( pageInd ) '.fig'] );
    %     print( '-dpsc', '-painters', '-append', '-cmyk', '-r300', [plotFilePathRoot '.ps'] );
    %     print( '-dpsc', '-painters', '-append', '-cmyk', '-r300', [plotFilePathRoot '_page' num2str( pageInd ) '.ps'] );
    %     close;
end

end
