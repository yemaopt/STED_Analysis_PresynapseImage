% Written by Ye Ma, Department of Biomedical Engineering, Johns Hopkins University

function Results = CalculateIntensityAndSize(image, points, threshold, windowSize, axPlot)

if nargin < 5
    axPlot = [];
end

[sizey, sizex] = size(image);
image = image - min(image(:));

Results = nan(size(points,1), 4);

if isempty(points)
    return
end

for spot_indx = 1:size(points,1)

    X0 = round(points(spot_indx,1));
    Y0 = round(points(spot_indx,2));

    xstart  = max(1, X0-windowSize);
    xfinish = min(sizex, X0+windowSize);
    ystart  = max(1, Y0-windowSize);
    yfinish = min(sizey, Y0+windowSize);

    img_tmp = double(image(ystart:yfinish, xstart:xfinish));

    localX = X0 - xstart + 1;
    localY = Y0 - ystart + 1;

    if localX < 1 || localX > size(img_tmp,2) || localY < 1 || localY > size(img_tmp,1)
        continue
    end

    peakVal = img_tmp(localY, localX);
    C = contourc(img_tmp, [threshold threshold] .* double(peakVal));
    S = contourdata(C);

    if isempty(S)
        continue
    end

    area_tmp = zeros(1, numel(S));
    for contour_indx = 1:numel(S)
        area_tmp(contour_indx) = polyarea(S(contour_indx).xdata, S(contour_indx).ydata);
    end

    [maxArea, id] = max(area_tmp);

    isClosed = hypot( ...
        S(id).xdata(1) - S(id).xdata(end), ...
        S(id).ydata(1) - S(id).ydata(end)) < 2;

    if ~isClosed
        continue
    end

    [xx, yy] = meshgrid(1:size(img_tmp,2), 1:size(img_tmp,1));
    indicator_in = inpolygon(xx(:), yy(:), S(id).xdata, S(id).ydata);
    indx_in = find(indicator_in);

    Results(spot_indx,1) = maxArea;
    Results(spot_indx,2) = sum(img_tmp(indx_in));
    Results(spot_indx,3) = X0;
    Results(spot_indx,4) = Y0;

    if ~isempty(axPlot) && isvalid(axPlot)
        hold(axPlot,'on')
        plot(axPlot, S(id).xdata + xstart - 1, ...
                     S(id).ydata + ystart - 1, ...
                     'LineWidth', 2)
        scatter(axPlot, xx(indx_in) + xstart - 1, ...
                        yy(indx_in) + ystart - 1, ...
                        8, 'filled')
        hold(axPlot,'off')
    end
end
end