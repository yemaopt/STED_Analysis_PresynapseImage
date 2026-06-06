function Results = CalculateIntensityAndSize(image, points, threshold, windowSize, axPlot)

if nargin < 5
    axPlot = [];
end

[sizey, sizex] = size(image);
image = double(image);
image = image - min(image(:));

% Results columns:
% 1 Area_px
% 2 TotalIntensity_bgCorrected
% 3 X0
% 4 Y0
% 5 EquivalentDiameter_px
% 6 PeakIntensity
% 7 LocalBackground
% 8 ThresholdValue
Results = nan(size(points,1), 8);

if isempty(points)
    return
end

for spot_indx = 1:size(points,1)

    X0 = round(points(spot_indx,1));
    Y0 = round(points(spot_indx,2));

    xstart  = max(1, X0 - windowSize);
    xfinish = min(sizex, X0 + windowSize);
    ystart  = max(1, Y0 - windowSize);
    yfinish = min(sizey, Y0 + windowSize);

    img_tmp = image(ystart:yfinish, xstart:xfinish);

    localX = X0 - xstart + 1;
    localY = Y0 - ystart + 1;

    if localX < 1 || localX > size(img_tmp,2) || ...
            localY < 1 || localY > size(img_tmp,1)
        continue
    end

    localBackground = median(img_tmp(:));
    peakIntensity = img_tmp(localY, localX);

    % Smooth lightly
    img_smooth = imgaussfilt(img_tmp, 0.7);

    % Use stricter threshold
    thresholdValue = localBackground + threshold * ...
        (peakIntensity - localBackground);

    mask = img_smooth >= thresholdValue;

    % Remove tiny noise
    mask = bwareaopen(mask, 3);

    % Keep only component containing selected peak
    CC = bwconncomp(mask, 8);

    if CC.NumObjects == 0
        continue
    end

    peakLinearIndex = sub2ind(size(img_tmp), localY, localX);

    componentID = [];
    for c = 1:CC.NumObjects
        if any(CC.PixelIdxList{c} == peakLinearIndex)
            componentID = c;
            break
        end
    end

    if isempty(componentID)
        continue
    end

    punctaIdx = CC.PixelIdxList{componentID};
    [yy, xx] = ind2sub(size(img_tmp), punctaIdx);
    distFromPeak = sqrt((xx - localX).^2 + (yy - localY).^2);
    maxRadius_px = 3.5;   % try 3–5 px for 30 nm pixels
    keepIdx = distFromPeak <= maxRadius_px;
    punctaIdx = punctaIdx(keepIdx);

    fprintf('Spot %d: %d pixels in puncta\n', ...
        spot_indx, ...
        numel(punctaIdx));

    Area_px = numel(punctaIdx);
    TotalIntensity = sum(img_tmp(punctaIdx) - localBackground);
    EquivalentDiameter_px = 2 * sqrt(Area_px / pi);

    Results(spot_indx,1) = Area_px;
    Results(spot_indx,2) = TotalIntensity;
    Results(spot_indx,3) = X0;
    Results(spot_indx,4) = Y0;
    Results(spot_indx,5) = EquivalentDiameter_px;
    Results(spot_indx,6) = peakIntensity;
    Results(spot_indx,7) = localBackground;
    Results(spot_indx,8) = thresholdValue;

    if ~isempty(axPlot) && isvalid(axPlot)

        punctaMask = false(size(img_tmp));
        punctaMask(punctaIdx) = true;

        B = bwboundaries(punctaMask,'noholes');

        hold(axPlot,'on')

        for b = 1:numel(B)
            boundaryY = B{b}(:,1) + ystart - 1;
            boundaryX = B{b}(:,2) + xstart - 1;

            plot(axPlot, boundaryX, boundaryY, ...
                'c-', ...
                'LineWidth', 3)
        end

        scatter(axPlot, X0, Y0, ...
            40, ...
            'mo', ...
            'LineWidth', 1.5)

        hold(axPlot,'off')
    end
end
end