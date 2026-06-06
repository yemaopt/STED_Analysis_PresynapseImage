% Written by Ye Ma, Department of Biomedical Engineering, Johns Hopkins University, 2022

clear
clc
close all
addpath ./functions
set(0,'DefaultFigureWindowStyle','docked')
selpath =  uigetdir('', 'Please select the folder containing the results');
imglist=dir(fullfile(selpath,sprintf('*.tif')));
load([selpath,'/','FinalResults.mat'], 'PunctaStatsResults', 'DistanceResults', 'ContourResults');

PixelSize = 30;
NumBin=8;
%%
Area=PunctaStatsResults(:,2).*PixelSize.*PixelSize./1e6;
DistanceToActiveZone=DistanceResults(:,1);
TotalIntensity=PunctaStatsResults(:,1);
PunctaParasList=table(DistanceToActiveZone,Area,TotalIntensity);
%%  Display the results ---- histogram
close all
figure;
histogram(DistanceToActiveZone,NumBin,'Normalization','probability');
xlabel('Distance from active zone boundary');
ylabel('Fraction of dynamin cluster')
set(gca,'FontSize',38);

figure;
histogram(Area,NumBin);
xlabel('Area of the puncta (nm^2)');
ylabel('Number of dynamin puncta')
set(gca,'FontSize',38);

figure;
histogram(TotalIntensity,NumBin);
xlabel('Total Intensity (A.U.)');
ylabel('Number of dynamin puncta')
set(gca,'FontSize',38);

% Display the results ---- scatter plot
figure;
scatter(DistanceToActiveZone,Area,100,'k','filled');
xlabel('Distance from the active zone (nm)');
ylabel('Area of the puncta (nm^2)')
set(gca,'FontSize',38);

figure;
scatter(DistanceToActiveZone,TotalIntensity,100,'k','filled');
xlabel('Distance from the active zone (nm)');
ylabel('Total Intensity (A.U.)')
set(gca,'FontSize',38);

%%
figure;
H=histogram(DistanceToActiveZone,NumBin,'Normalization','probability');
BinEdges=H.BinEdges';
%BinEdges=[-80   -30    20    70   120   170   220   270   320]';

H=histogram(DistanceToActiveZone,BinEdges,'Normalization','probability');
xlabel('Distance from active zone boundary');
ylabel('Fraction of dynamin cluster')
yyaxis left

Mean_Area=zeros(length(BinEdges)-1,1);
Mean_TotalIntensity=zeros(length(BinEdges)-1,1);
SE_Area=zeros(length(BinEdges)-1,1);
SE_TotalIntensity=zeros(length(BinEdges)-1,1);
X=zeros(length(BinEdges)-1,1);

puncta_index_set=cell(length(BinEdges)-1,1);

for ii=1:length(BinEdges)-1
    ind = find(DistanceToActiveZone>=BinEdges(ii) & DistanceToActiveZone<BinEdges(ii+1));
    Mean_Area(ii) = nanmean(Area(ind));
    Mean_TotalIntensity(ii) = nanmean(TotalIntensity(ind));
    SE_Area(ii) = nanstd(Area(ind));
    SE_TotalIntensity(ii) = nanstd(TotalIntensity(ind));
    X(ii)=(BinEdges(ii)+BinEdges(ii+1))/2;
    puncta_index_set{ii}=ind;
end
hold on
yyaxis right
p2=errorbar(X,Mean_Area,SE_Area,'c-','LineWidth',2);
set(gca,'FontSize',30);
ylabel('Puncta size (nm)')
legend([p2],'Area')
%%
figure;
histogram(DistanceToActiveZone,BinEdges,'Normalization','probability');
xlabel('Distance from active zone boundary');
ylabel('Fraction of dynamin cluster')
yyaxis left
hold on
yyaxis right
p4=errorbar(X,Mean_TotalIntensity,SE_TotalIntensity,'r-','LineWidth',2);
set(gca,'FontSize',30);
ylabel('Total Intensity (A.U.)')
legend([p4],'Intensity')
%%
BinEdgeStart=BinEdges(1:end-1);
BinEdgeEnd=BinEdges(2:end);
Fraction=H.Values';
BinStatsList=table(BinEdgeStart,BinEdgeEnd,Fraction,Mean_Area,SE_Area,Mean_TotalIntensity,SE_TotalIntensity,puncta_index_set);
%%
BinPunctaAreaList=zeros(500,length(BinEdges)-1);
BinPunctaIntensityList=zeros(500,length(BinEdges)-1);
BinPunctaDistanceList=zeros(500,length(BinEdges)-1);

for ii=1:length(BinEdges)-1
    for jj=1:length(puncta_index_set{ii})
        BinPunctaAreaList(jj,ii)=PunctaParasList.Area(puncta_index_set{ii}(jj));
        BinPunctaIntensityList(jj,ii)=PunctaParasList.TotalIntensity(puncta_index_set{ii}(jj));
        BinPunctaDistanceList(jj,ii)=PunctaParasList.DistanceToActiveZone(puncta_index_set{ii}(jj));
    end
end
