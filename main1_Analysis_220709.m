% Written by Ye Ma, Department of Biomedical Engineering, Johns Hopkins University, 2022

clear
clc
close all
addpath ./functions
set(0,'DefaultFigureWindowStyle','docked')
selpath = uigetdir('', 'Please select the folder containing the ROI images');
imglist=dir(fullfile(selpath,sprintf('*.tif')));
mkdir([selpath,'\screenshot'])
mkdir([selpath,'\output'])

% For dynamin channel
PeakThresholdFactor=0.35;
V=0.5:0.05:0.9;
PixelSize=30;
windowSize=5;
preview_flag=0;

distance_accumulated=[];
NumActiveZone=0;
NumDynCluster=0;
PunctaStatsResults=[];
DistanceResults=[];
ContourResults={};

%% Selection
for imgind=1:length(imglist)
    Tmp = loadtiff([imglist(imgind).folder,'/',imglist(imgind).name]);
    Ch0_dyn=Tmp(:,:,1);
    Ch1_bas=Tmp(:,:,2);
    %% Preview
    [PeakY,PeakX,PeakVal]=localMaximum(double(Ch0_dyn),3,1,PeakThresholdFactor.*max(Ch0_dyn(:)));
    figure(1);
    annotation('textbox', [0.35, 0.98, 0, 0], 'string', imglist(imgind).name, 'Interpreter', 'none');
    subplot(341);hold on;imagesc(Ch0_dyn);scatter(PeakX,PeakY,'c*');axis equal;axis off;hold off;title('Ch dynamin preview');set(gca,'XDir','normal');set(gca,'YDir','reverse');
    C=contourc(double(Ch1_bas),V.*double(max(Ch1_bas(:))));
    subplot(342);imagesc(Ch1_bas);axis equal;axis off;hold on;title('Ch bassoon preview');set(gca,'XDir','normal');set(gca,'YDir','reverse');
    S=contourdata(C);
    for k=1:length(S)
        plot(S(k).xdata,S(k).ydata,'LineWidth',2)
    end
    hold off
    colormap(hot)
    subplot(343);imshowpair(Ch0_dyn,Ch1_bas);title('Composite preview');set(gca,'XDir','normal');set(gca,'YDir','reverse');
    %% Select dynamin cluster
    subplot(345);hold on;imagesc(Ch0_dyn);scatter(PeakX,PeakY,'c*');axis equal;axis off;title('Select dynamin cluster, end with ENTER');set(gca,'XDir','normal');set(gca,'YDir','reverse');
    [xi,yi]=getpts();
    if isempty(xi)
        close figure 1
        continue
    end
    NumDynCluster=NumDynCluster+length(xi);
    xlist=[];
    ylist=[];
    PeakList=PeakX+1j*PeakY;
    for i=1:length(xi)
        temp=abs(PeakList-xi(i)-1j.*yi(i));
        [~,ind]=min(temp);
        xlist(i)=PeakX(ind);
        ylist(i)=PeakY(ind);
    end
    xlist=xlist';
    ylist=ylist';
    scatter(xlist,ylist,'ko')
    hold off;
    %% Select active zone contour
    subplot(346);imagesc(Ch1_bas);axis equal;axis off;hold on;title('Select active zone boundary, end with ENTER');set(gca,'XDir','normal');set(gca,'YDir','reverse');
    S=contourdata(C);
    for k=1:length(S)
        plot(S(k).xdata,S(k).ydata,'LineWidth',2)
    end
    colormap(hot)
    [x_SelectedContour,y_SelectedContour]=getpts();
    if isempty(x_SelectedContour)
        close figure 1
        continue
    end
    NumActiveZone=NumActiveZone+length(x_SelectedContour);
    distance_current_set=[];xy_set=[];ind_selectedContour_set=[];
    for jjj = 1:length(x_SelectedContour)
        d=100.*ones(length(S),1);
        for k=1:length(S)
            [~,distance_tmp,~] = distance2curve([S(k).xdata,S(k).ydata],[x_SelectedContour(jjj),y_SelectedContour(jjj)],'linear');
            d(k)=distance_tmp;
        end
        [~,ind_selectedContour]=min(d);
        ind_selectedContour_set(jjj)=ind_selectedContour;
        %%  Calculate distance
        [xy,distance_current_tmp,t_a] = distance2curve([S(ind_selectedContour).xdata,S(ind_selectedContour).ydata],[xlist,ylist],'spline');
        distance_current_set(:,jjj)=distance_current_tmp;
        xy_set(:,:,jjj)=xy;
    end
    [distance_current,indmin]=min(distance_current_set,[],2);
    xy_current=[];
    for kkkk=1:length(indmin)
        xy_current(kkkk,:)=xy_set(kkkk,:,indmin(kkkk));
    end
    subplot(347);imshowpair(Ch0_dyn,Ch1_bas);axis equal;axis off;hold on;set(gca,'XDir','normal');set(gca,'YDir','reverse');
    for kkk=1:length(ind_selectedContour_set)
        plot(S(ind_selectedContour_set(kkk)).xdata,S(ind_selectedContour_set(kkk)).ydata,'c-','LineWidth',2)
    end
    scatter(xlist,ylist,'c*')
    scatter(xy_current(:,1),xy_current(:,2),'r*')
    line([xlist,xy_current(:,1)]',[ylist,xy_current(:,2)]','color',[0 0 1])
    title('Distance calculated')
    hold off
    ActiveZoneContour={};
    %%  Display the distance result
    for iii=1:length(distance_current)
        if inpolygon(xlist(iii),ylist(iii),S(ind_selectedContour_set(indmin(iii))).xdata,S(ind_selectedContour_set(indmin(iii))).ydata)
           distance_current(iii)=-distance_current(iii);
        end
        ActiveZoneContour{iii}=S(ind_selectedContour_set(indmin(iii)));
    end
    distance_current=distance_current.*PixelSize;
    subplot(348);stem(distance_current);title('Distance distribution');xlabel('Dynamin cluster #');ylabel('Distance from active zone boundary')
    distance_accumulated=[distance_accumulated',distance_current']';
    subplot(344);histogram(distance_accumulated,10,'Normalization','probability');title('Accumulated distance distribution');xlabel('Distance from active zone boundary');ylabel('Fraction of dynamin cluster')
    str1={'# of active zone = ',num2str(NumActiveZone)};
    str2={'# of dynamin cluster = ',num2str(NumDynCluster)};
    annotation('textbox', [0.91, 0.91, 0.8, 0], 'string', str1);
    annotation('textbox', [0.91, 0.85, 0.8, 0], 'string', str2);
    text(2,7,str1)
    text(2,10,str2)
    %% Gaussian fitting for dynamin size
    points = [xlist,ylist];
   
    figure(1);subplot(3,4,9);imagesc(Ch0_dyn);colormap(hot);axis equal;axis off
    Results_temp = CalculateIntensityAndSize(Ch0_dyn,points,0.75,windowSize);
  % Results_temp = Results_temp(~isnan(Results_temp(:,1)),:);
    npts=size(Results_temp,1);
    
    for i2=1:npts
        text(Results_temp(i2,3),Results_temp(i2,4),num2str(i2),'Color','cyan','FontSize',14)
        text(Results_temp(i2,3),Results_temp(i2,4)+2,num2str(Results_temp(i2,1)),'Color','cyan','FontSize',14)
        text(Results_temp(i2,3),Results_temp(i2,4)+4,num2str(Results_temp(i2,2)),'Color','cyan','FontSize',14)
    end
    
    Results_temp(:,5)=imgind;
    distance_current(:,2)=imgind;
    PunctaStatsResults=[PunctaStatsResults',Results_temp']';
    DistanceResults=[DistanceResults',distance_current']';
    ContourResults=[ContourResults,ActiveZoneContour];
    Area=PunctaStatsResults(1,:);
    TotalIntensity=PunctaStatsResults(2,:);
    %%
    subplot(3,4,10);histogram(PunctaStatsResults(:,1).*PixelSize.*PixelSize./1e6,10);title('Histogram of Area');xlabel('Area (um^2)');ylabel('Count');
    subplot(3,4,11);histogram(PunctaStatsResults(:,2),10);title('Histogram of Intensity');xlabel('Total Intensity (A.U.)');ylabel('Count');
    pause()
    %%  Save
    %screenshot
    print('-dtiff','-r100',[imglist(imgind).name,'_screenshot.tiff']);
    movefile([imglist(imgind).name,'_screenshot.tiff'], [selpath,'\screenshot'])
    %save
    save([imglist(imgind).name,'_result.mat'], 'Results_temp', 'distance_current', 'ActiveZoneContour');
    movefile([imglist(imgind).name,'_result.mat'], [selpath,'\output'])
    %%  Close figure
    close figure 1
end
save([selpath,'\FinalResults.mat'], 'PunctaStatsResults', 'DistanceResults', 'ContourResults','Area','TotalIntensity');

