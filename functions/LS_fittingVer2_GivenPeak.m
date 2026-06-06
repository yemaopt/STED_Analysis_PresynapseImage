% Output
% Results = [Amp,xPos,yPos,sigmaX,sigmaY,background=0,X0,Y0,sum(ccd_image(:))]
%%
function [Results] = LS_fittingVer2_GivenPeak(image,points,threshold,windowSize,preview_flag)
epsilon1 =  10^-5;
epsilon2 =  10^-5;
maxIter  = 500;
curvefitoptions = optimset( 'lsqcurvefit');
curvefitoptions = optimset( curvefitoptions,'Jacobian' ,'off','Display', 'off',  'TolX', epsilon2, 'TolFun', epsilon1,'MaxPCGIter',1,'MaxIter',maxIter);
Results= GaussFit(image,points,threshold,curvefitoptions,windowSize,preview_flag);
end

function fitresults = GaussFit(image,points,threshold,curvefitoptions,windowSize,preview_flag) % do the storm calculation

[sizey,sizex] = size(image);
image=image-min(image(:));

% [PeakY,PeakX,~] = localMaximum(image,3,1,threshold.*max(image(:)));
% points=[PeakX,PeakY];

if isempty(points)==0
    nrspots=length(points(:,1));
    fitresults=zeros(nrspots,10);
    
    for i=1:nrspots
        X0=points(i,1);
        Y0=points(i,2);
        
        xstart =  X0-windowSize;
        xfinish = X0+windowSize;
        ystart =  Y0-windowSize;
        yfinish = Y0+windowSize;
        
        xstart=max(1,xstart);ystart=max(1,ystart);
        xfinish=min(sizex,xfinish);yfinish=min(sizey,yfinish);
        
        X_POSim = X0-xstart+1;
        Y_POSim = Y0-ystart+1;
        
        img=double(image(ystart:yfinish,xstart:xfinish));

        background = 0;
        brightness = max(img(:));
        widthStart = 1.5;
        
        initguess=double([brightness,X_POSim,Y_POSim,widthStart,widthStart,background]);
        xLim = [X_POSim-0.001 X_POSim+0.001];
        yLim = [Y_POSim-0.001 Y_POSim+0.001];
        sigmaLim = [0.5  6];
        
        [fitParams,res] = Gauss2d_Fit(img,initguess,xLim,yLim,sigmaLim,curvefitoptions);
        
        % assign the data
        amplitude = fitParams(1);
        xPos = fitParams(2)+xstart-1;
        yPos = fitParams(3)+ystart-1;
        sigmaX = max(fitParams(4),fitParams(5));
        sigmaY = min(fitParams(4),fitParams(5));
        rotAngle = fitParams(6);
        background = fitParams(7);
        fitresults(i,:)=[amplitude,xPos,yPos,sigmaX,sigmaY,rotAngle,background,X0,Y0,sum(img(:))];
    end
    
    if preview_flag==1
        for i=1:nrspots
            X0=points(i,1);
            Y0=points(i,2);
            xstart =  X0-windowSize;
            xfinish = X0+windowSize;
            ystart =  Y0-windowSize;
            yfinish = Y0+windowSize;
            xstart=max(1,xstart);ystart=max(1,ystart);
            xfinish=min(sizex,xfinish);yfinish=min(sizey,yfinish);
            img_tmp=double(image(ystart:yfinish,xstart:xfinish));
            figure(101)
            subplot(3,3,i);imagesc(img_tmp);
        end
    end
    
else
    fitresults=[];
end
end

%----------------------------------------------------------------------------------------------

function [fitParam,res] = Gauss2d_Fit(inputIm,initguess,xLim,yLim,sigmaLim,curvefitoptions)
Astart = initguess(1);
xStart = initguess(2);
yStart = initguess(3);
widthStartX = initguess(4);
widthStartY = initguess(5);
BGstart = initguess(6);

xMin = xLim(1);
xMax = xLim(2);
yMin = yLim(1);
yMax = yLim(2);

sigmaMin = sigmaLim(1);
sigmaMax = sigmaLim(2);
[sizey,sizex] = size(inputIm);
[X,Y]= meshgrid(1:sizex,1:sizey);    % combine X and Y into a single matrix, and split them in Gauss2d
grid = [X Y];

initGuess5Vector = [Astart xStart yStart widthStartX widthStartY 0 BGstart];

lb = [1       xMin    yMin    sigmaMin    sigmaMin  -pi/4   0     ];
ub = [65535   xMax    yMax    sigmaMax    sigmaMax   pi/4   65535 ];

try
    [fitParam, res] = ...
        lsqcurvefit(@(x, xdata) Gauss2d(x, xdata), ...
        initGuess5Vector ,grid ,inputIm ,...
        lb,ub,curvefitoptions);
catch ME
    if strcmp(ME.identifier,'optim:snls:InvalidUserFunction')
        fitParam = [0 0 0 0 0 0 0];
        res = 0;
    else
        rethrow(ME);
    end
end

if fitParam(1)< 0
    fitParam = [0 0 0 0 0 0 0];
end
end


%%
function F = Gauss2d(a, data)
% a(1) - A
% a(2) - Xpos
% a(3) - Ypos
% a(4) - sigmaX
% a(5) - sigmaY
% a(6) - Rotation Angle
% a(7) - B
% a(6) = 0;
a(7) = 0;

x0rot = a(2)*cos(a(6)) - a(3)*sin(a(6));
y0rot = a(2)*sin(a(6)) + a(3)*cos(a(6));

[~,sizex] = size(data);
sizex= sizex/2;
X = data(:,1:sizex);
Y = data(:,sizex+1:end);

Xrot = X*cos(a(6)) - Y*sin(a(6));
Yrot = X*sin(a(6)) + Y*cos(a(6));

expPart = exp(-(Xrot-x0rot).^2./(2.*(a(4).^2))-(Yrot-y0rot).^2./(2.*(a(5).^2)));
F =  a(1)*expPart + a(7);
end
