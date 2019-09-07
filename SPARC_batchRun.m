if 0 
if 0
    
    if 0
        clear all
        close all
        
        titley = 'FlpOut'
        bouties = [1]
        SPARC_batch_FlpOut
        SPARC_batchSetup
        
        clear all
        close all
        
        titley = 'SPARC'
        bouties = [1]
        SPARC_batch_SPARC
        SPARC_batchSetup
    end
    
    %%
    
    clear all
    close all
    
    titley = 'FlpOut'
    bouties = [1:3]
    colorz = 'm'
    SPARC_batch_FlpOut
    SPARC_batchSetup
    
    clearvars -except mapsShifted_cum_mean_one mapsShifted_cum_std_half_one tuningInd_cum_one
    close all
    
    titley = 'SPARC'
    bouties = [1:3]
    colorz = 'k'
    SPARC_batch_SPARC
    SPARC_batchSetup
    
    close all
    
    figure
    errorbar(squeeze(mapsShifted_cum_mean_one), mapsShifted_cum_std_half_one, 'm', 'LineWidth', 2);
    hold on
    errorbar(squeeze(mapsShifted_cum_mean), mapsShifted_cum_std_half, 'k', 'LineWidth', 2);
    title({[titley ' bouts: ' num2str(bouties)]; ['mean tuning curve' ' DSI: ' num2str(tcParams_new.tuningInd)]})
    ylim([0 1.3])
    xlim([0 9])
    hold off
    
    figure
    scatter(ones(1, length(tuningInd_cum_one)) .* 0.6, [tuningInd_cum_one], 'm');
    hold on
    scatter(ones(1, length(tuningInd_cum)).* 0.2, [tuningInd_cum], 'k');
    scatter(0.6, (mean(tuningInd_cum_one)), 'filled', 'm');
    scatter(0.2, (mean(tuningInd_cum)), 'filled', 'k');
    title({[titley ' bouts: ' num2str(bouties)]; 'DSI distribution'})
    ylim([0 1])
    xlim([0 2])
    hold off
    
    SaveOpenFigures ('/scratch/users/wienecke/figures/SPARC_figure_190612/');
    
    
    %%
    
    clear all
    close all
    
    titley = 'FlpOut'
    bouties = [1:5]
    colorz = 'm'
    SPARC_batch_FlpOut
    SPARC_batchSetup
    
    clearvars -except mapsShifted_cum_mean_one mapsShifted_cum_std_half_one tuningInd_cum_one
    close all
    
    titley = 'SPARC'
    bouties = [1:5]
    colorz = 'k'
    SPARC_batch_SPARC
    SPARC_batchSetup
    
    close all
    
    figure
    errorbar(squeeze(mapsShifted_cum_mean_one), mapsShifted_cum_std_half_one, 'm', 'LineWidth', 2);
    hold on
    errorbar(squeeze(mapsShifted_cum_mean), mapsShifted_cum_std_half, 'k', 'LineWidth', 2);
    title({[titley ' bouts: ' num2str(bouties)]; ['mean tuning curve' ' DSI: ' num2str(tcParams_new.tuningInd)]})
    ylim([0 1.3])
    xlim([0 9])
    hold off
    
    figure
    scatter(ones(1, length(tuningInd_cum_one)) .* 0.6, [tuningInd_cum_one], 'm');
    hold on
    scatter(ones(1, length(tuningInd_cum)).* 0.2, [tuningInd_cum], 'k');
    scatter(0.6, (mean(tuningInd_cum_one)), 'filled', 'm');
    scatter(0.2, (mean(tuningInd_cum)), 'filled', 'k');
    title({[titley ' bouts: ' num2str(bouties)]; 'DSI distribution'})
    ylim([0 1])
    xlim([0 2])
    hold off
    
    SaveOpenFigures ('/scratch/users/wienecke/figures/SPARC_figure_190612/');
    
    
    %%
    
    clear all
    close all
    
    titley = 'FlpOut_short'
    bouties = [1:3]
    colorz = 'm'
    SPARC_batch_FlpOut_short
    SPARC_batchSetup
    
    clearvars -except mapsShifted_cum_mean_one mapsShifted_cum_std_half_one tuningInd_cum_one
    close all
    
    titley = 'SPARC'
    bouties = [1:3]
    colorz = 'k'
    SPARC_batch_SPARC
    SPARC_batchSetup
    
    close all
    
    figure
    errorbar(squeeze(mapsShifted_cum_mean_one), mapsShifted_cum_std_half_one, 'm', 'LineWidth', 2);
    hold on
    errorbar(squeeze(mapsShifted_cum_mean), mapsShifted_cum_std_half, 'k', 'LineWidth', 2);
    title({[titley ' bouts: ' num2str(bouties)]; ['mean tuning curve short' ' DSI: ' num2str(tcParams_new.tuningInd)]})
    ylim([0 1.3])
    xlim([0 9])
    hold off
    
    figure
    scatter(ones(1, length(tuningInd_cum_one)) .* 0.6, [tuningInd_cum_one], 'm');
    hold on
    scatter(ones(1, length(tuningInd_cum)).* 0.2, [tuningInd_cum], 'k');
    scatter(0.6, (mean(tuningInd_cum_one)), 'filled', 'm');
    scatter(0.2, (mean(tuningInd_cum)), 'filled', 'k');
    title({[titley ' bouts: ' num2str(bouties)]; 'DSI distribution short'})
    ylim([0 1])
    xlim([0 2])
    hold off
    
    SaveOpenFigures ('/scratch/users/wienecke/figures/SPARC_figure_190612/');
    
    
    %%
    
end

clear all
close all

singleFigures = false
saveFigures = true
local = false

titley = 'FlpOut_short'
bouties = [1:5]
colorz = 'k'
SPARC_batch_FlpOut_short
SPARC_batchSetup

clearvars -except singleFigures saveFigures local mapsShifted_cum_mean_one mapsShifted_cum_std_half_one tuningInd_cum_one
close all

titley = 'SPARC'
bouties = [1:5]
colorz = 'm'
SPARC_batch_SPARC
SPARC_batchSetup
%%

widthy=3
figure
errorbar(squeeze(mapsShifted_cum_mean_one), mapsShifted_cum_std_half_one, 'k', 'LineWidth', widthy);
hold on
errorbar(squeeze(mapsShifted_cum_mean), mapsShifted_cum_std_half, 'm', 'LineWidth', widthy);
title({[titley ' bouts: ' num2str(bouties)]; ['mean tuning curve short' ' DSI: ' num2str(tcParams_new.tuningInd)]})
ylim([0 1.3])
xlim([0 9])
hold off

saveDir = ['/scratch/users/wienecke/SPARC_data/'];

if ~isdir( saveDir )
    mkdir( saveDir )
end
save( [saveDir 'mapsShifted_cum_mean_one_'], 'mapsShifted_cum_mean_one', '-v7.3', '-mat' );
save( [saveDir 'mapsShifted_cum_mean_'], 'mapsShifted_cum_mean', '-v7.3', '-mat' );

%%

[h,p] = ttest2(tuningInd_cum,tuningInd_cum_one,'Alpha', 0.01,'Vartype','unequal')
sizey=60
figure
scatter(ones(1, length(tuningInd_cum_one)) .* 0.8, [tuningInd_cum_one], sizey, 'k');
hold on
scatter(ones(1, length(tuningInd_cum)).* 0.4, [tuningInd_cum], sizey, 'm');
scatter(0.8, (mean(tuningInd_cum_one)),sizey, 'filled', 'k');
scatter(0.4, (mean(tuningInd_cum)), sizey, 'filled', 'm');
title({[titley ' bouts: ' num2str(bouties) ' h/p: ' num2str(h) '_' num2str(p)]; 'DSI distribution short'})
ylim([0 1])
xlim([0 6])
hold off

saveDir = ['/scratch/users/wienecke/SPARC_data/'];

if ~isdir( saveDir )
    mkdir( saveDir )
end
save( [saveDir 'tuningInd_cum_one_'], 'tuningInd_cum_one', '-v7.3', '-mat' );
save( [saveDir 'tuningInd_cum_'], 'tuningInd_cum', '-v7.3', '-mat' );


%%

if saveFigures
    SaveOpenFigures ('/scratch/users/wienecke/figures/SPARC_figure_190612/');
end

clearvars -except singleFigures saveFigures local
clc;
end
local = false 
saveFigures = true 
singleFigures = false 

%set up common variables
functionalChannel = 1;
plotFlag = 0;
generateFitData = false;
upSampleFactor = 1;
labelStruct = [];
roiStruct = [];

if local
    basePath = '~/Documents/ambrose/leprechaunMat/181019.0.cfrw';
    flyInd = 1;
    roiInd = [1];
else
    basePath = '/scratch/users/wienecke/leprechaunMat/190428.0.cfrw';
    flyInd = 9;
    roiInd = 2;
end

[imagingStruct, registrationStruct, maskStruct, timingStruct, matchStruct, stimParamsMetadata, stimParams, tuningStruct, roiInfo] = ...
    LoadSession(basePath, flyInd, 1, true);
imagingStruct.stackRaw = imagingStruct.stackRaw(:, roiInd, :);

[taImagingStruct, taTimingStruct] = TrialAverage( imagingStruct, registrationStruct, timingStruct, matchStruct, functionalChannel, 1, [], [] );
[~, ~, epochStacks1] = PlotRoiTimeCourses_magenta( taImagingStruct, taTimingStruct, [], [], 1, 'dF_F0', [], [], '','','','','','','','', 2, [34/255 160/255 31/255]);

saveDir = ['/scratch/users/wienecke/SPARC_data/'];

if ~isdir( saveDir )
    mkdir( saveDir )
end

save( [saveDir 'epochStacks1_'], 'epochStacks1', '-v7.3', '-mat' );

if saveFigures
    SaveOpenFigures ('/scratch/users/wienecke/figures/SPARC_figure_190612/');
end


clearvars -except singleFigures saveFigures local
clc;

%set up common variables
functionalChannel = 1;
plotFlag = 0;
generateFitData = false;
upSampleFactor = 1;
labelStruct = [];
roiStruct = [];

if local
    basePath = '~/Documents/ambrose/leprechaunMat/181019.0.cfrw';
    flyInd = 1;
    roiInd = [1];
else
    basePath = '/scratch/users/wienecke/leprechaunMat/190415.0.cfrw';
    flyInd = 0;
    roiInd = 9;
end
[imagingStruct, registrationStruct, maskStruct, timingStruct, matchStruct, stimParamsMetadata, stimParams, tuningStruct, roiInfo] = ...
    LoadSession(basePath, flyInd, 1, true);
imagingStruct.stackRaw = imagingStruct.stackRaw(:, roiInd, :);

[taImagingStruct, taTimingStruct] = TrialAverage( imagingStruct, registrationStruct, timingStruct, matchStruct, functionalChannel, 1, [], [] );
[~, ~, epochStacks2] = PlotRoiTimeCourses_magenta( taImagingStruct, taTimingStruct, [], [], 1, 'dF_F0', [], [], '','','','','','','','', 2, [0 0 0]);
saveDir = ['/scratch/users/wienecke/SPARC_data/'];

if ~isdir( saveDir )
    mkdir( saveDir )
end

save( [saveDir 'epochStacks2_'], 'epochStacks2', '-v7.3', '-mat' );

if saveFigures
    SaveOpenFigures ('/scratch/users/wienecke/figures/SPARC_figure_190612/');
end


%%

a1=[50.30 35.94 43.00 44.07 40.55 43.88]
a2=[16.28 12.26 9.21 10.42 13.66 8.75]
a3=[6.56 3.50 8.77 3.53 2.57 2.82]

figure
scatter(ones(1, length(a1)) .* 0.5, [a1], sizey, 'c');
hold on
scatter(ones(1, length(a2)).* 1, [a2], sizey, 'b');
scatter(ones(1, length(a3)).* 1.5, [a3], sizey, 'm');

scatter(0.5, (mean(a1)),sizey, 'filled', 'c');
scatter(1, (mean(a2)),sizey, 'filled', 'b');

scatter(1.5, (mean(a3)), sizey, 'filled', 'm');
title(['SPARC counts---means: ' num2str(mean(a1)) ' ' num2str(mean(a2)) ' ' num2str(mean(a3))])
ylim([0 100])
xlim([0 6])
hold off

%% 

clear all 
close all

sizey = 60 

one = [50.0 50.7 38.7 52.2 43.9]
two = [44.74 44.62 53.97 48.53 46.43]
three = [7.69 21.54 15.07 18.31 19.70]
four = [25.00 22.73 18.82 17.54 16.83]
five = [16.67 9.23 12.00 19.35 15.38]
six = [2.67 3.07 2.93 4.93 3.33]

figure
scatter(ones(1, length(one)) .* 0.5, [one], sizey, 'c');
hold on
scatter(ones(1, length(two)).* 1, [two], sizey, 'b');
scatter(ones(1, length(three)).* 1.5, [three], sizey, 'm');
scatter(ones(1, length(four)).* 2, [four], sizey, 'r');
scatter(ones(1, length(five)).* 2.5, [five], sizey, 'g');
scatter(ones(1, length(six)).* 3, [six], sizey, 'k');

scatter(0.5, (mean(one)),sizey, 'filled', 'c');
scatter(1, (mean(two)),sizey, 'filled', 'b');
scatter(1.5, (mean(three)),sizey, 'filled', 'm');
scatter(2.0, (mean(four)),sizey, 'filled', 'r');
scatter(2.5, (mean(five)),sizey, 'filled', 'g');
scatter(3.0, (mean(six)),sizey, 'filled', 'k');

title(['SPARC CD8GFP counts---means: ' num2str(mean(one)) ' ' num2str(mean(two)) ' ' num2str(mean(three)) ' ' num2str(mean(four)) ' ' num2str(mean(five)) ' ' num2str(mean(six))])
ylim([0 100])
xlim([0 6])
hold off

    SaveOpenFigures ('~/Documents');
