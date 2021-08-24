%FAKE WAVE PLAN VIEW AND CROSS SECTION PLOT
savePath = 'D:\Dropbox\matlab\RadarSim\CrossSection\4panel\';

%load conus

%MAKE WAVE FIELDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
volumeTime = 60*4;%in seconds

wave1Length = 30720;%in meters

wave1Perturbation = 4;%in m/s

wave1Dir = 270;

wave1Height = 4500;

wave1Speed = 8;%in m/s

justPlotting = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if ~justPlotting
[modelX,modelY,modelZ] = meshgrid(linspace(-74.9,-70.8,600),linspace(39.2,42.5,600),linspace(0,8000,200));
%[xEast,yNorth,zUp]=geodetic2enu(modelY,modelX,modelZ,22.5,-32.5,0,referenceEllipsoid('earth', 'm'));

%create fake radar
[ radarStruct ] = createRadarStructure(-72.859, 40.861, 1, 0.5, 200000, 250, 35.9);%KOKX

[Xdist,Xaz] = distance(40.861,-72.859,40.861,radarStruct.lonRadar,referenceEllipsoid('earth','km'));Xdist(Xaz>180)=Xdist(Xaz>180).*-1;
[Ydist,Yaz] = distance(40.861,-72.859,radarStruct.latRadar,-72.859,referenceEllipsoid('earth','km'));Ydist(Yaz>90)=Ydist(Yaz>90).*-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=0:15
    
    %Make wave schematics
    [Uwave1,Vwave1,Wwave1] = waveGenerator(modelX,modelY,modelZ,wave1Length,wave1Height,wave1Dir,wave1Perturbation,ii*volumeTime*wave1Speed,0);
    [Uwave2,Vwave2,Wwave2] = waveGenerator(modelX,modelY,modelZ,wave1Length,wave1Height,wave1Dir,wave1Perturbation,ii*volumeTime*wave1Speed,85.236);%4min at 8 m/s
    [ interpolantsWave1 ] = createInterpolants( modelX,modelY,modelZ,Uwave1,Vwave1,Wwave1);
    [ DVwave1 ] = calculateDV( radarStruct, interpolantsWave1);
    [ interpolantsWave2 ] = createInterpolants( modelX,modelY,modelZ,Uwave2,Vwave2,Wwave2);
    [ DVwave2 ] = calculateDV( radarStruct, interpolantsWave2);
    
    %end
    
    %PLOTTING
    
    domainWidth = distance(40.85,-74.9,40.85,-70.8,referenceEllipsoid('earth','km'));
    
    FH=figure('Position',[10 10 1920 840]);
    
    AH1 = subplot(4,2,1);
    %AH1 = axes('Position',[1/16 6.5/8 5/16 1/8]);
    imagesc([-1.*domainWidth./2 domainWidth./2],[0 8000],squeeze(Wwave1(300,:,:))');
    set(AH1,'YDir','normal','FontSize',14);
    ylabel('Height (m)','FontSize',14);
    xlabel('Distance From Radar (km)','FontSize',14);
    ylim([0 6000]);xlim([0 50]);
    colormap(AH1,circshift(fliplr(Colormap_DV_RdBu(64)),2,2));
    caxis([-0.235 0.235]);colorbar;
    title('Idealized Wave Simulated Vertical Velocity (m/s)','FontSize',14,'Interpreter','none');
    
    AH2 = subplot(4,2,2);
    %AH2 = axes('Position',[10/16 6.5/8 5/16 1/8]);
    imagesc([-1.*domainWidth./2 domainWidth./2],[0 8000],squeeze(Wwave2(300,:,:))');
    set(AH2,'YDir','normal','FontSize',14);
    ylabel('Height (m)','FontSize',14);
    xlabel('Distance From Radar (km)','FontSize',14);
    ylim([0 6000]);xlim([0 50]);
    colormap(AH2,circshift(fliplr(Colormap_DV_RdBu(64)),2,2));
    caxis([-0.235 0.235]);colorbar;
    title('Idealized Wave Simulated Vertical Velocity (m/s)','FontSize',14,'Interpreter','none');
    
    AH3 = subplot(4,2,[3 5 7]);
    %AH3 = axes('Position',[1/16 0.5/8 5/16 5/8]);
    surf(squeeze(Xdist),squeeze(Ydist),squeeze(DVwave1),'EdgeColor','none','FaceColor','texturemap');view(2);axis image;
    set(AH3,'box','on','XGrid','off','YGrid','off','ZGrid','off','FontSize',14,'Layer','top','xtick',-150:50:150,'ytick',-150:50:150);
    ylabel('Distance (km)','FontSize',14); xlabel('Distance (km)','FontSize',14);
    colormap(AH3,flipud(Colormap_DV_RdBu(64)));
    caxis([-5 5]);colorbar;
    hold on
    plot3([0 50],[0 0],[10 10],'k','LineWidth',2);
    %title('Idealized Wave Simulated Radial Velocity (m/s)','FontSize',14,'Interpreter','none');
    
    AH4 = subplot(4,2,[4 6 8]);
    %AH4 = axes('Position',[10/16 0.5/8 5/16 5/8]);
    surf(squeeze(Xdist),squeeze(Ydist),squeeze(DVwave2),'EdgeColor','none','FaceColor','texturemap');view(2);axis image;
    set(AH4,'box','on','XGrid','off','YGrid','off','ZGrid','off','FontSize',14,'Layer','top','xtick',-150:50:150,'ytick',-150:50:150);
    ylabel('Distance (km)','FontSize',14); xlabel('Distance (km)','FontSize',14);
    colormap(AH4,flipud(Colormap_DV_RdBu(64)));
    caxis([-5 5]);colorbar;
    hold on
    plot3([0 50],[0 0],[10 10],'k','LineWidth',2);
    %title('Idealized Wave Simulated Radial Velocity (m/s)','FontSize',14,'Interpreter','none');
    
    FH.PaperUnits = 'inches';
    FH.PaperPosition = [0 0 16 7];
    print([savePath 'WaveTilt_xsect_4panel_frame' num2str(ii,'%02u') '.png'],'-dpng','-r120');
    close
    
end




writerObj = VideoWriter([savePath 'WaveTilt_xsect_animation'],'MPEG-4');
writerObj.FrameRate = 6;
open(writerObj);
frameList=dir([savePath 'WaveTilt_xsect_4panel_frame*.png']);
for ii=1:length(frameList)
    Img=imread([savePath frameList(ii).name]);
    writeVideo(writerObj,Img);
end
close(writerObj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Horizontal Velocity - Time 1
% FH3=figure('Position',[10 10 1920 1080]);
% imagesc([-1.*domainWidth./2 domainWidth./2],[0 8000],squeeze(Uwave1(300,:,:))');
% set(gca,'YDir','normal','FontSize',16);
% ylabel('Height (m)','FontSize',16);
% xlabel('Distance From Radar (km)','FontSize',16);
% ylim([0 6000]);xlim([0 50]);
% colormap(gca,Colormap_DV_RdBu(64));
% caxis([-5 5]);colorbar;
% title('Idealized Wave Horizontal Velocity (m/s) - Time 0:00 minutes','FontSize',20,'Interpreter','none');
%
% FH3.PaperUnits = 'inches';
% FH3.PaperPosition = [0 0 16 4.5];
% %print([savePath 'IdealizedWave_xsect_U_time1.png'],'-dpng','-r120');
% %close
%
% %Horizontal Velocity - Time 2
% FH4=figure('Position',[10 10 1920 1080]);
% imagesc([-1.*domainWidth./2 domainWidth./2],[0 8000],squeeze(Uwave2(300,:,:))');
% set(gca,'YDir','normal','FontSize',16);
% ylabel('Height (m)','FontSize',16);
% xlabel('Distance From Radar (km)','FontSize',16);
% ylim([0 6000]);xlim([0 50]);
% colormap(gca,Colormap_DV_RdBu(64));
% caxis([-5 5]);colorbar;
% title('Idealized Wave Horizontal Velocity (m/s) - Time +4:00 minutes','FontSize',20,'Interpreter','none');
%
% FH4.PaperUnits = 'inches';
% FH4.PaperPosition = [0 0 16 4.5];
% %print([savePath 'IdealizedWave_xsect_U_time2.png'],'-dpng','-r120');
% %close
%
% %Instantaneous Horizontal Divergence - Time 1
% FH5=figure('Position',[10 10 1920 1080]);
% imagesc([-1.*domainWidth./2 domainWidth./2],[0 8000],squeeze(circshift(Uwave1(300,:,:),1,2) - circshift(Uwave1(300,:,:),-1,2))');
% set(gca,'YDir','normal','FontSize',16);
% ylabel('Height (m)','FontSize',16);
% xlabel('Distance From Radar (km)','FontSize',16);
% ylim([0 6000]);xlim([0 50]);
% colormap(gca,rot90(Colormap_DV_RdBu(64),2));
% caxis([-2 2]);colorbar;
% title('Idealized Wave Horizontal Convergence (m/s) - Time 0:00 minutes','FontSize',20,'Interpreter','none');
%
% FH5.PaperUnits = 'inches';
% FH5.PaperPosition = [0 0 16 4.5];
% %print([savePath 'IdealizedWave_xsect_HorzDiv_time1.png'],'-dpng','-r120');
% %close
%
% %Instantaneous Horizontal Divergence - Time 2
% FH6=figure('Position',[10 10 1920 1080]);
% imagesc([-1.*domainWidth./2 domainWidth./2],[0 8000],squeeze(circshift(Uwave2(300,:,:),1,2) - circshift(Uwave2(300,:,:),-1,2))');
% set(gca,'YDir','normal','FontSize',16);
% ylabel('Height (m)','FontSize',16);
% xlabel('Distance From Radar (km)','FontSize',16);
% ylim([0 6000]);xlim([0 50]);
% colormap(gca,rot90(Colormap_DV_RdBu(64),2));
% caxis([-2 2]);colorbar;
% title('Idealized Wave Horizontal Convergence (m/s) - Time +4:00 minutes','FontSize',20,'Interpreter','none');
%
% FH6.PaperUnits = 'inches';
% FH6.PaperPosition = [0 0 16 4.5];
% %print([savePath 'IdealizedWave_xsect_HorzDiv_time2.png'],'-dpng','-r120');
%
% %close
%
% %Horizontal Diff
% FH7=figure('Position',[10 10 1920 1080]);
% imagesc([-1.*domainWidth./2 domainWidth./2],[0 8000],squeeze(Uwave2(300,:,:))' - squeeze(Uwave1(300,:,:))');
% set(gca,'YDir','normal','FontSize',16);
% ylabel('Height (m)','FontSize',16);
% xlabel('Distance From Radar (km)','FontSize',16);
% ylim([0 6000]);xlim([0 50]);
% colormap(gca,Colormap_DV_RdBu(64));
% caxis([-5 5]);colorbar;
% title('Idealized Wave Horizontal Velocity Difference (m/s)','FontSize',20,'Interpreter','none');
%
% FH7.PaperUnits = 'inches';
% FH7.PaperPosition = [0 0 16 4.5];
% %print([savePath 'IdealizedWave_xsect_HorzDifference_time1.png'],'-dpng','-r120');
% %close
%
% %Pseudo Wave ID - pure horizontal case
% FH8=figure('Position',[10 10 1920 1080]);
% imagesc([-1.*domainWidth./2 domainWidth./2],[0 8000],squeeze(Uwave2(300,:,:))' - squeeze(Uwave1(300,:,:))');
% set(gca,'YDir','normal','FontSize',16);
% ylabel('Height (m)','FontSize',16);
% xlabel('Distance From Radar (km)','FontSize',16);
% ylim([0 6000]);xlim([0 50]);
% colormap(gca,[0 0 0; 1 1 1]);
% caxis([-0.065 -0.060]);colorbar;
% title('Idealized Wave ID','FontSize',20,'Interpreter','none');
%
% FH8.PaperUnits = 'inches';
% FH8.PaperPosition = [0 0 16 4.5];
% %print([savePath 'IdealizedWave_xsect_WaveID.png'],'-dpng','-r120');
% %close
%
%
%
% FH11=figure;
% surf(squeeze(Xdist),squeeze(Ydist),squeeze(DVwave2) - squeeze(DVwave1),'EdgeColor','none','FaceColor','texturemap');view(2);axis image;
% set(gca,'box','on','XGrid','off','YGrid','off','ZGrid','off','FontSize',16);
% colormap(gca,[0 0 0; 1 1 1]);
% caxis([-0.065 -0.060]);colorbar;
% title('Idealized Wave Simulated Radial Velocity Wave ID','FontSize',20,'Interpreter','none');
%
% FH11.PaperUnits = 'inches';
% FH11.PaperPosition = [0 0 9 9];
% %print([savePath 'IdealizedWave_planview_WaveID.png'],'-dpng','-r120');
% %close
%
% FH12=figure;
% surf(squeeze(Xdist),squeeze(Ydist),squeeze(DVwave2) - squeeze(DVwave1),'EdgeColor','none','FaceColor','texturemap');view(2);axis image;
% set(gca,'box','on','XGrid','off','YGrid','off','ZGrid','off','FontSize',16);
% colormap(gca,flipud(Colormap_DV_RdBu(64)));
% caxis([-5 5]);colorbar;
% title('Idealized Wave Simulated Radial Velocity Difference (m/s)','FontSize',20,'Interpreter','none');
%
% FH12.PaperUnits = 'inches';
% FH12.PaperPosition = [0 0 9 9];
% %print([savePath 'IdealizedWave_planview_DVdiff.png'],'-dpng','-r120');
% %close

%
% FH=figure;
% imagesc([min(squeeze(modelX(300,:,1))) max(squeeze(modelX(300,:,1)))],[min(squeeze(modelZ(300,3,:))) max(squeeze(modelZ(300,3,:)))],squeeze(Uwave1(300,:,:))');
% set(gca,'YDir','normal');
% ylim([0 6000]);xlim([-73 -72.5]);
% colormap(gca,imadjust(Colormap_DV_BWR(64,0.25,0.75),[0 1],[0.7 1]));colorbar;
% hold on;
% [C1,h1] = contour(squeeze(modelX(300,:,:))',squeeze(modelZ(300,:,:))',squeeze(Wwave1(300,:,:))',0.02:0.04:0.26,'-k');
% [C2,h2] = contour(squeeze(modelX(300,:,:))',squeeze(modelZ(300,:,:))',squeeze(Wwave1(300,:,:))',-0.26:0.04:0.02,':k');
% clabel(C1,h1,0.02:0.08:0.26);
% clabel(C2,h2,-0.26:0.08:0.02);
% xlabel('Longitude','FontSize',16);
% ylabel('Height (m)','FontSize',16);
% title('Idealized Wave Cross Section','FontSize',20,'Interpreter','none');
% FH.PaperUnits = 'inches';
% FH.PaperPosition = [0 0 16 9];
%print([savePath 'IdealizedWave_crosssection.png'],'-dpng','-r120');

%close




