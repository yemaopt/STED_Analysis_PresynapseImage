% Output
% Results = [Amp,xPos,yPos,sigmaX,sigmaY,background=0,X0,Y0,sum(ccd_image(:))]
%%
function [Results] = CalculateIntensityAndSize(image,points,threshold,windowSize)

[sizey,sizex] = size(image);
image=image-min(image(:));

if isempty(points)==0
    nrspots=length(points(:,1));
    
    for spot_indx=1:nrspots
        X0=points(spot_indx,1);
        Y0=points(spot_indx,2);
        xstart =  X0-windowSize;
        xfinish = X0+windowSize;
        ystart =  Y0-windowSize;
        yfinish = Y0+windowSize;
        xstart=max(1,xstart);ystart=max(1,ystart);
        xfinish=min(sizex,xfinish);yfinish=min(sizey,yfinish);
        img_tmp=double(image(ystart:yfinish,xstart:xfinish));
        
        C=contourc(double(img_tmp),[threshold,threshold].*double(img_tmp(X0-xstart+1,Y0-ystart+1)));
        S=contourdata(C);
        
        area_tmp=zeros(1,length(S));
        for contour_indx=1:length(S)
            area_tmp(contour_indx)=polyarea(S(contour_indx).xdata,S(contour_indx).ydata);
        end
        
        [~,id]=max(area_tmp);
        
        if ((S(id).xdata(1)-S(id).xdata(end))^2+(S(id).ydata(1)-S(id).ydata(end))^2)^0.5<2
            
            [yy,xx]=meshgrid(1:size(img_tmp,2),1:size(img_tmp,1));
            indicator_in=inpolygon(yy(:),xx(:),S(id).xdata,S(id).ydata);
            indx_in=find(indicator_in==1);
            [indx_row,indx_col] = ind2sub(size(img_tmp,1),indx_in);
            
            Results(spot_indx,1)=max(area_tmp);
            Results(spot_indx,2)=sum(img_tmp(indx_in));
            Results(spot_indx,3)=(xstart+xfinish)/2;
            Results(spot_indx,4)=(ystart+yfinish)/2;
            
            figure(1)
            subplot(349);hold on;axis equal;axis off;colormap(hot)
            plot(S(id).xdata+xstart-1,S(id).ydata+ystart-1,'LineWidth',2)
            scatter(indx_col+xstart-1,indx_row+ystart-1)
            hold off
            
%             figure(101)
%             subplot(3,3,spot_indx);imagesc(img_tmp);hold on
%             plot(S(id).xdata,S(id).ydata,'LineWidth',2)
%             scatter(indx_col,indx_row)
%             hold off
        end
        
        if ((S(id).xdata(1)-S(id).xdata(end))^2+(S(id).ydata(1)-S(id).ydata(end))^2)^0.5>2
            
            Results(spot_indx,1)=NaN;
            Results(spot_indx,2)=NaN;
            Results(spot_indx,3)=NaN;
            Results(spot_indx,4)=NaN;
 
        end
        
    end
end
end

