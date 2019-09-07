function [imagingStruct, registrationStruct, maskStruct, timingStruct, matchStruct, stimParamsMetadata, stimParams, tuningStruct, roiInfo] = ...
    LoadSession(basePath, flyInd, stackInd, useImagingRoiStruct, roiStruct_uniqueID)
% LOADSESSION - Load all data structures from a particular day for one stack of one fly.
%
% For example, to load data for imaging session 1 for fly0 in 140515.0.jcsl you would run:
% imagingStruct, registrationStruct, timingStruct, matchStruct, stimParamsMetadata, stimParams] = ...
%     LoadSession('/hsgs/projects/ganguli/leprechaunMat/140515.0.jcsl/', 0, 1)
%
% stackInd can be either a numerical index or it can be a string to search
% for. 170204 leongjcs

if ischar( stackInd )
    stackInd = GetStackInd( basePath, flyInd, stackInd );
end

if ~exist('useImagingRoiStruct', 'var') || isempty(useImagingRoiStruct)
    useImagingRoiStruct = false;
end

if ~exist( 'roiStruct_uniqueID', 'var' ) || isempty( roiStruct_uniqueID )
    roiStruct_uniqueID = '';
end

baseName = FindSessionPrefix(basePath, flyInd);

if useImagingRoiStruct
    try
        imagingStruct = LoadStruct(basePath, flyInd, stackInd, 'imagingRoiStruct', roiStruct_uniqueID);
        roiInfo = LoadStruct(basePath, flyInd, stackInd, 'imagingRoiStruct_roiInfo', roiStruct_uniqueID);
    catch
      error( 'No imagingRoiStruct with specified roiStruct uniqueID.');
    end
else
    try
        imagingStructFilename = [baseName 'imagingStruct.mat'];
        f = matfile(imagingStructFilename);
        imagingStruct = f.imagingStruct(1, stackInd);
        clear f;
    catch
        disp( 'WARNING: No imagingStruct found. imagingStruct returned is an empty array.' );
        
        imagingStruct = [];
    end
    roiInfo = [];
end

try
    registrationStructFilename = [baseName 'registrationStruct.mat'];
    f = matfile(registrationStructFilename);
    registrationStruct = f.registrationStruct(1, stackInd);
    clear f;
catch
    disp( 'WARNING: No registrationStruct found. registrationStruct returned is an empty array.' );
    
    registrationStruct = [];
end

try
    maskStructFilename = [baseName 'maskStruct.mat'];
    f = matfile(maskStructFilename);
    maskStruct = f.maskStruct;
    clear f;
catch
    disp( 'WARNING: No maskStruct found. maskStruct returned is an empty array.' );
    
    maskStruct = [];
end

try
    matchStructFilename = [baseName 'matchStruct.mat'];
    f = matfile(matchStructFilename);
    matchStruct = f.matchStruct;
    clear f;
catch
    disp( 'WARNING: No matchStruct found. matchStruct returned is an empty array.' );
    
    matchStruct = [];
end

try
    timingStructFilename = [baseName 'timingStruct.mat'];
    f = matfile(timingStructFilename);
    if matchStruct.timingStructInds(stackInd) ~= 0
        timingStruct = f.timingStruct(1, matchStruct.timingStructInds(stackInd));
    else
        timingStruct = f.timingStruct;
    end
    clear f;
catch
    disp( 'WARNING: No timingStruct found. timingStruct returned is an empty array.' );
    
    timingStruct = [];
end

try
    stimulusStructsFilename = [baseName 'stimulusStructs.mat'];
    f = matfile(stimulusStructsFilename);
    if matchStruct.timingStructInds(stackInd) ~= 0
        stimParams = f.stimParams(1, matchStruct.timingStructInds(stackInd));
        stimParams = stimParams{1};
        stimParamsMetadata = f.stimParamsMetadata(1, matchStruct.timingStructInds(stackInd));
    else
        stimParams = f.stimParams;
        stimParamsMetadata = f.stimParamsMetadata;
    end
    clear f;
catch
    disp( 'WARNING: No stimulusStructs found. stimulusStructs returned are empty arrays.' );
    
    stimParams = [];
    stimParamsMetadata = [];
end

try
    tuningStructFilename = [baseName 'tuningStruct.mat'];
    f = matfile(tuningStructFilename);
    tuningStruct = f.tuningStruct;
    clear f;
catch
    disp( 'WARNING: No tuningStruct found. tuningStruct returned is an empty array.' );
    
    tuningStruct = [];
end

try
    if matchStruct.timingStructInds(stackInd) ~= 0
        matchStruct = IndexIntoMatchStruct( matchStruct, stackInd );
    else
        disp( 'WARNING: imagingStruct did not match with other structs. Loaded full matchStruct and all timingStructs, stimParams, and stimParamsMetadata.' );
    end
catch
end

end
