% Written by Ye Ma, Department of Biomedical Engineering, Johns Hopkins University
% Optimized interface/structure/readability version.
%
% Notes:
% - Requires helper functions in ./functions:
%   loadtiff, localMaximum, contourdata, distance2curve, CalculateIntensityAndSize
% - Uses tiledlayout/nexttile, available in newer MATLAB versions.

clear
clc
close all
warning('off', 'images:initSize:adjustingMag');

%% -------------------- Configuration --------------------
cfg = defaultConfig();

addpath(cfg.functionsPath)

selpath = uigetdir('', 'Please select the folder containing the ROI images');
if isequal(selpath,0)
    error('No folder selected. Analysis cancelled.');
end

imglist = dir(fullfile(selpath,'*.tif'));
if isempty(imglist)
    error('No .tif files found in the selected folder: %s', selpath);
end

outputDirs = prepareOutputFolders(selpath);

%% -------------------- Initialize results --------------------
results = initializeResults();

%% -------------------- Main image loop --------------------
for imgind = 1:numel(imglist)

    imgInfo = imglist(imgind);
    imgPath = fullfile(imgInfo.folder, imgInfo.name);

    fprintf('\nProcessing image %d/%d: %s\n', imgind, numel(imglist), imgInfo.name);

    try
        imgData = loadImageChannels(imgPath);
    catch ME
        warning('Failed to load image %s. Skipping. Error: %s', imgInfo.name, ME.message);
        continue
    end

    Ch0_dyn = imgData.Ch0_dyn;
    Ch1_bas = imgData.Ch1_bas;

    [PeakX, PeakY] = detectDynaminPeaks(Ch0_dyn, cfg);
    contourStruct = detectBassoonContours(Ch1_bas, cfg);

    if isempty(PeakX)
        warning('No dynamin peaks detected in %s. Skipping.', imgInfo.name);
        continue
    end

    if isempty(contourStruct)
        warning('No bassoon contours detected in %s. Skipping.', imgInfo.name);
        continue
    end

    fig = createAnalysisFigure(imgInfo.name, cfg);
    t = tiledlayout(fig, 3, 4, ...
        'TileSpacing', cfg.tileSpacing, ...
        'Padding', cfg.tilePadding);

    plotPreviewPanels(t, Ch0_dyn, Ch1_bas, PeakX, PeakY, contourStruct, imgInfo.name);

    %% -------------------- Manual selection: dynamin clusters --------------------
    axSelectDyn = nexttile(t,5);
    plotDynaminSelectionPanel(axSelectDyn, Ch0_dyn, PeakX, PeakY);

    [xi, yi] = getpts(axSelectDyn);
    if isempty(xi)
        fprintf('No dynamin clusters selected. Skipping %s.\n', imgInfo.name);
        close(fig)
        continue
    end

    [xlist, ylist] = mapSelectedPointsToNearestPeaks(xi, yi, PeakX, PeakY);
    hold(axSelectDyn,'on')
    scatter(axSelectDyn, xlist, ylist, cfg.selectedMarkerSize, 'ko', 'LineWidth', 1.5)
    hold(axSelectDyn,'off')

    %% -------------------- Manual selection: active-zone contours --------------------
    axSelectAZ = nexttile(t,6);
    plotContourSelectionPanel(axSelectAZ, Ch1_bas, contourStruct);

    [x_SelectedContour, y_SelectedContour] = getpts(axSelectAZ);
    if isempty(x_SelectedContour)
        fprintf('No active-zone contour selected. Skipping %s.\n', imgInfo.name);
        close(fig)
        continue
    end

    %% -------------------- Distance calculation --------------------
    [distance_current, xy_current, activeZoneContour, selectedContourIdx] = ...
        calculateClusterToContourDistances(xlist, ylist, x_SelectedContour, y_SelectedContour, contourStruct, cfg);

    if isempty(distance_current)
        warning('Distance calculation returned empty result for %s. Skipping.', imgInfo.name);
        close(fig)
        continue
    end

    results.NumActiveZone = results.NumActiveZone + numel(x_SelectedContour);
    results.NumDynCluster = results.NumDynCluster + numel(xlist);

    axDistance = nexttile(t,7);
    plotDistancePanel(axDistance, Ch0_dyn, Ch1_bas, contourStruct, selectedContourIdx, ...
        xlist, ylist, xy_current);

    axStem = nexttile(t,8);
    plotCurrentDistanceDistribution(axStem, distance_current);

    results.distance_accumulated = [results.distance_accumulated; distance_current(:)];

    axAccum = nexttile(t,4);
    plotAccumulatedDistanceDistribution(axAccum, results.distance_accumulated);

    axSummary = nexttile(t,12);
    addSummaryText(axSummary, results.NumActiveZone, results.NumDynCluster);

    %% -------------------- Puncta size/intensity calculation --------------------
    points = [xlist(:), ylist(:)];

    axPuncta = nexttile(t,9);
    imagesc(axPuncta, Ch0_dyn)
    formatImageAxis(axPuncta)
    colormap(axPuncta, hot)
    title(axPuncta, 'Dynamin puncta size/intensity')

    Results_temp = CalculateIntensityAndSize(Ch0_dyn, points, cfg.gaussianFitThreshold, cfg.windowSize);
    Results_temp = annotatePunctaResults(axPuncta, Results_temp);

    Results_temp(:,5) = imgind;
    distance_with_img_index = [distance_current(:), imgind * ones(numel(distance_current),1)];

    results.PunctaStatsResults = [results.PunctaStatsResults; Results_temp];
    results.DistanceResults = [results.DistanceResults; distance_with_img_index];
    results.ContourResults = [results.ContourResults, activeZoneContour];

    axArea = nexttile(t,10);
    plotAreaHistogram(axArea, results.PunctaStatsResults, cfg);

    axIntensity = nexttile(t,11);
    plotIntensityHistogram(axIntensity, results.PunctaStatsResults);

    drawnow

    fprintf('Review the figure. Press any key in MATLAB command window to continue.\n');
    pause()

    %% -------------------- Save per-image outputs --------------------
    saveImageOutputs(fig, imgInfo.name, outputDirs, Results_temp, distance_with_img_index, activeZoneContour, cfg);

    close(fig)
end

%% -------------------- Save final outputs --------------------
if isempty(results.PunctaStatsResults)
    Area = [];
    TotalIntensity = [];
else
    Area = results.PunctaStatsResults(:,1);
    TotalIntensity = results.PunctaStatsResults(:,2);
end

PunctaStatsResults = results.PunctaStatsResults;
DistanceResults = results.DistanceResults;
ContourResults = results.ContourResults;

save(fullfile(selpath,'FinalResults.mat'), ...
    'PunctaStatsResults', 'DistanceResults', 'ContourResults', ...
    'Area', 'TotalIntensity');

fprintf('\nAnalysis complete.\nFinal results saved to:\n%s\n', fullfile(selpath,'FinalResults.mat'));

%% ========================================================================
% Local helper functions
% ========================================================================

function cfg = defaultConfig()
    cfg.functionsPath = './functions';

    cfg.PeakThresholdFactor = 0.35;
    cfg.PixelSize_nm = 30;
    cfg.windowSize = 5;
    cfg.gaussianFitThreshold = 0.75;

    cfg.contourLevels = 0.5:0.05:0.9;

    cfg.figurePosition = [80 80 1550 950];
    cfg.figureColor = 'w';
    cfg.tileSpacing = 'compact';
    cfg.tilePadding = 'compact';
    cfg.selectedMarkerSize = 50;

    cfg.screenshotResolution = '-r150';
end

function outputDirs = prepareOutputFolders(selpath)
    outputDirs.screenshot = fullfile(selpath,'screenshot');
    outputDirs.output = fullfile(selpath,'output');

    if ~exist(outputDirs.screenshot,'dir')
        mkdir(outputDirs.screenshot)
    end

    if ~exist(outputDirs.output,'dir')
        mkdir(outputDirs.output)
    end
end

function results = initializeResults()
    results.distance_accumulated = [];
    results.NumActiveZone = 0;
    results.NumDynCluster = 0;
    results.PunctaStatsResults = [];
    results.DistanceResults = [];
    results.ContourResults = {};
end

function imgData = loadImageChannels(imgPath)
    Tmp = loadtiff(imgPath);

    if ndims(Tmp) < 3 || size(Tmp,3) < 2
        error('Input image must contain at least two channels.');
    end

    imgData.Ch0_dyn = Tmp(:,:,1);
    imgData.Ch1_bas = Tmp(:,:,2);
end

function [PeakX, PeakY] = detectDynaminPeaks(Ch0_dyn, cfg)
    threshold = cfg.PeakThresholdFactor * double(max(Ch0_dyn(:)));
    [PeakY, PeakX, ~] = localMaximum(double(Ch0_dyn), 3, 1, threshold);
end

function contourStruct = detectBassoonContours(Ch1_bas, cfg)
    C = contourc(double(Ch1_bas), cfg.contourLevels .* double(max(Ch1_bas(:))));
    contourStruct = contourdata(C);
end

function fig = createAnalysisFigure(imageName, cfg)
    fig = figure( ...
        'Name', ['STED analysis: ', imageName], ...
        'NumberTitle', 'off', ...
        'Units', 'pixels', ...
        'Position', cfg.figurePosition, ...
        'Color', cfg.figureColor, ...
        'MenuBar', 'figure', ...
        'ToolBar', 'figure');
end

function plotPreviewPanels(t, Ch0_dyn, Ch1_bas, PeakX, PeakY, contourStruct, imageName)
    ax1 = nexttile(t,1);
    imagesc(ax1, Ch0_dyn)
    hold(ax1,'on')
    scatter(ax1, PeakX, PeakY, 25, 'c*')
    hold(ax1,'off')
    formatImageAxis(ax1)
    colormap(ax1, hot)
    title(ax1, 'Dynamin preview')

    ax2 = nexttile(t,2);
    imagesc(ax2, Ch1_bas)
    hold(ax2,'on')
    plotContours(ax2, contourStruct, 1.5)
    hold(ax2,'off')
    formatImageAxis(ax2)
    colormap(ax2, hot)
    title(ax2, 'Bassoon contours')

    ax3 = nexttile(t,3);
    imshowpair(Ch0_dyn, Ch1_bas, 'Parent', ax3)
    formatImageAxis(ax3)
    title(ax3, 'Composite preview')

    sgtitle(t, imageName, 'Interpreter', 'none', 'FontWeight', 'bold')
end

function plotDynaminSelectionPanel(ax, Ch0_dyn, PeakX, PeakY)
    imagesc(ax, Ch0_dyn)
    hold(ax,'on')
    scatter(ax, PeakX, PeakY, 25, 'c*')
    hold(ax,'off')
    formatImageAxis(ax)
    colormap(ax, hot)
    title(ax, {'Select dynamin clusters', 'Press ENTER when finished'})
end

function plotContourSelectionPanel(ax, Ch1_bas, contourStruct)
    imagesc(ax, Ch1_bas)
    hold(ax,'on')
    plotContours(ax, contourStruct, 1.5)
    hold(ax,'off')
    formatImageAxis(ax)
    colormap(ax, hot)
    title(ax, {'Select active-zone boundaries', 'Press ENTER when finished'})
end

function [xlist, ylist] = mapSelectedPointsToNearestPeaks(xi, yi, PeakX, PeakY)
    PeakList = PeakX(:) + 1i * PeakY(:);
    xlist = zeros(numel(xi),1);
    ylist = zeros(numel(yi),1);

    for i = 1:numel(xi)
        distanceToPeaks = abs(PeakList - xi(i) - 1i * yi(i));
        [~, ind] = min(distanceToPeaks);
        xlist(i) = PeakX(ind);
        ylist(i) = PeakY(ind);
    end
end

function [distance_nm, xy_current, activeZoneContour, selectedContourIdx] = ...
    calculateClusterToContourDistances(xlist, ylist, xSelected, ySelected, contourStruct, cfg)

    nSelectedContours = numel(xSelected);
    nClusters = numel(xlist);

    distance_current_set = nan(nClusters, nSelectedContours);
    xy_set = nan(nClusters, 2, nSelectedContours);
    selectedContourIdx = zeros(1, nSelectedContours);

    for j = 1:nSelectedContours
        d = inf(numel(contourStruct),1);

        for k = 1:numel(contourStruct)
            curveXY = [contourStruct(k).xdata, contourStruct(k).ydata];
            [~, distance_tmp, ~] = distance2curve(curveXY, [xSelected(j), ySelected(j)], 'linear');
            d(k) = distance_tmp;
        end

        [~, contourIdx] = min(d);
        selectedContourIdx(j) = contourIdx;

        curveXY = [contourStruct(contourIdx).xdata, contourStruct(contourIdx).ydata];
        [xy, distance_tmp, ~] = distance2curve(curveXY, [xlist(:), ylist(:)], 'spline');

        distance_current_set(:,j) = distance_tmp(:);
        xy_set(:,:,j) = xy;
    end

    [distance_pixels, minIdx] = min(distance_current_set, [], 2);

    xy_current = nan(nClusters,2);
    activeZoneContour = cell(1,nClusters);

    for i = 1:nClusters
        xy_current(i,:) = xy_set(i,:,minIdx(i));

        nearestContourIdx = selectedContourIdx(minIdx(i));
        activeZoneContour{i} = contourStruct(nearestContourIdx);

        insideContour = inpolygon( ...
            xlist(i), ylist(i), ...
            contourStruct(nearestContourIdx).xdata, ...
            contourStruct(nearestContourIdx).ydata);

        if insideContour
            distance_pixels(i) = -distance_pixels(i);
        end
    end

    distance_nm = distance_pixels(:) * cfg.PixelSize_nm;
end

function plotDistancePanel(ax, Ch0_dyn, Ch1_bas, contourStruct, selectedContourIdx, xlist, ylist, xy_current)
    imshowpair(Ch0_dyn, Ch1_bas, 'Parent', ax)
    hold(ax,'on')

    selectedContourIdx = unique(selectedContourIdx);
    for k = 1:numel(selectedContourIdx)
        idx = selectedContourIdx(k);
        plot(ax, contourStruct(idx).xdata, contourStruct(idx).ydata, 'c-', 'LineWidth', 2)
    end

    scatter(ax, xlist, ylist, 35, 'c*')
    scatter(ax, xy_current(:,1), xy_current(:,2), 35, 'r*')

    for i = 1:numel(xlist)
        line(ax, [xlist(i), xy_current(i,1)], [ylist(i), xy_current(i,2)], ...
            'Color', [0 0 1], 'LineWidth', 1)
    end

    hold(ax,'off')
    formatImageAxis(ax)
    title(ax, 'Distance calculated')
end

function plotCurrentDistanceDistribution(ax, distance_current)
    stem(ax, distance_current, 'filled')
    title(ax, 'Current image distances')
    xlabel(ax, 'Dynamin cluster #')
    ylabel(ax, 'Distance from boundary (nm)')
    grid(ax, 'on')
end

function plotAccumulatedDistanceDistribution(ax, distance_accumulated)
    histogram(ax, distance_accumulated, 10, 'Normalization', 'probability')
    title(ax, 'Accumulated distance distribution')
    xlabel(ax, 'Distance from boundary (nm)')
    ylabel(ax, 'Fraction of clusters')
    grid(ax, 'on')
end

function Results_temp = annotatePunctaResults(ax, Results_temp)
    if isempty(Results_temp)
        return
    end

    for i = 1:size(Results_temp,1)
        text(ax, Results_temp(i,3), Results_temp(i,4), num2str(i), ...
            'Color', 'cyan', 'FontSize', 12, 'FontWeight', 'bold')

        text(ax, Results_temp(i,3), Results_temp(i,4)+2, sprintf('A=%.2f', Results_temp(i,1)), ...
            'Color', 'cyan', 'FontSize', 10)

        text(ax, Results_temp(i,3), Results_temp(i,4)+4, sprintf('I=%.2f', Results_temp(i,2)), ...
            'Color', 'cyan', 'FontSize', 10)
    end
end

function plotAreaHistogram(ax, PunctaStatsResults, cfg)
    area_um2 = PunctaStatsResults(:,1) .* cfg.PixelSize_nm .* cfg.PixelSize_nm ./ 1e6;
    histogram(ax, area_um2, 10)
    title(ax, 'Area distribution')
    xlabel(ax, 'Area (um^2)')
    ylabel(ax, 'Count')
    grid(ax, 'on')
end

function plotIntensityHistogram(ax, PunctaStatsResults)
    histogram(ax, PunctaStatsResults(:,2), 10)
    title(ax, 'Intensity distribution')
    xlabel(ax, 'Total intensity (A.U.)')
    ylabel(ax, 'Count')
    grid(ax, 'on')
end

function addSummaryText(ax, numActiveZone, numDynCluster)
    cla(ax)
    axis(ax,'off')

    text(ax, 0.05, 0.75, 'Summary', ...
        'Units','normalized', ...
        'FontSize',14, ...
        'FontWeight','bold')

    text(ax, 0.05, 0.55, sprintf('# active zones: %d', numActiveZone), ...
        'Units','normalized', ...
        'FontSize',12)

    text(ax, 0.05, 0.40, sprintf('# dynamin clusters: %d', numDynCluster), ...
        'Units','normalized', ...
        'FontSize',12)
end

function saveImageOutputs(fig, imageName, outputDirs, Results_temp, distance_current, ActiveZoneContour, cfg)
    screenshotName = [imageName, '_screenshot.tiff'];
    screenshotPath = fullfile(outputDirs.screenshot, screenshotName);
    print(fig, '-dtiff', cfg.screenshotResolution, screenshotPath)

    resultName = [imageName, '_result.mat'];
    resultPath = fullfile(outputDirs.output, resultName);
    save(resultPath, 'Results_temp', 'distance_current', 'ActiveZoneContour');
end

function plotContours(ax, contourStruct, lineWidth)
    for k = 1:numel(contourStruct)
        plot(ax, contourStruct(k).xdata, contourStruct(k).ydata, 'LineWidth', lineWidth)
    end
end

function formatImageAxis(ax)
    axis(ax, 'image')
    axis(ax, 'off')
    set(ax, 'XDir', 'normal')
    set(ax, 'YDir', 'reverse')
end
