%Two-wave, no-background demo

%savePath = 'D:\Dropbox\matlab\RadarSim\Wave_XSect_5panel';

%MAKE WAVE FIELDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveFigures = false;

justPlotting = false;

plotVector = false;

volumeTime = 60*4;%in seconds

wave1Length = 25000;%in meters

wave1Perturbation = 4;%in m/s

wave1Dir = 290;%deg from North

wave1Height = 4500;%meters

wave1Speed = 15;%in m/s

wave2Length = 8000;%in meters

wave2Perturbation = 2;%in m/s

wave2Dir = 250;%deg from North

wave2Height = 1000;%meters

wave2Speed = 8;%in m/s

radarRange = 150000;%metres

rangePadding = 5000;%meters

domainHeight = 7500;%meters

radarLoc = [-72.859, 40.861]; % lon, lat % KOKX

xSectDist = 40;%km
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~justPlotting
    % compute wave grid extents
    [waveEdgeN,~] = reckon(radarLoc(2),radarLoc(1),radarRange + rangePadding,0,referenceEllipsoid('earth', 'm'));
    [~,waveEdgeE] = reckon(radarLoc(2),radarLoc(1),radarRange + rangePadding,90,referenceEllipsoid('earth', 'm'));
    [waveEdgeS,~] = reckon(radarLoc(2),radarLoc(1),radarRange + rangePadding,180,referenceEllipsoid('earth', 'm'));
    [~,waveEdgeW] = reckon(radarLoc(2),radarLoc(1),radarRange + rangePadding,270,referenceEllipsoid('earth', 'm'));
    
        
    [modelX,modelY,modelZ] = meshgrid(linspace(waveEdgeW,waveEdgeE,400),linspace(waveEdgeS,waveEdgeN,400),linspace(0,domainHeight,150));
    %[xEast,yNorth,zUp]=geodetic2enu(modelY,modelX,modelZ,22.5,-32.5,0,referenceEllipsoid('earth', 'm'));
    
    %create fake radar
    [ radarStruct ] = createRadarStructure(radarLoc(1), radarLoc(2), 1, 0.5, radarRange, 1000, 35.9);%KOKX
    
    %Radar Parameters
    Tilts=[0.5 0.9 1.3 1.8 2.4 3.1 4.0 5.1 6.4 8.0 10.0 12.0 14.0 16.7 19.5];%VCP12
    [ radarStruct_2 ] = createRadarStructure(radarLoc(1), radarLoc(2), 15, Tilts, radarRange, 1000, 35.9);%KOKX
    
    [Xdist,Xaz] = distance(radarLoc(2),radarLoc(1),radarLoc(2),radarStruct.lonRadar,referenceEllipsoid('earth','km'));Xdist(Xaz>180)=Xdist(Xaz>180).*-1;
    [Ydist,Yaz] = distance(radarLoc(2),radarLoc(1),radarStruct.latRadar,radarLoc(1),referenceEllipsoid('earth','km'));Ydist(Yaz>90)=Ydist(Yaz>90).*-1;
    
    [Xdist2,Xaz2] = distance(radarLoc(2),radarLoc(1),radarLoc(2),radarStruct_2.lonRadar,referenceEllipsoid('earth','km'));Xdist2(Xaz2>180)=Xdist2(Xaz2>180).*-1;
    [Ydist2,Yaz2] = distance(radarLoc(2),radarLoc(1),radarStruct_2.latRadar,radarLoc(1),referenceEllipsoid('earth','km'));Ydist2(Yaz2>90)=Ydist2(Yaz2>90).*-1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=0:29
    
    [Uwave1,Vwave1,Wwave1] = waveGenerator(modelX,modelY,modelZ,wave1Length,wave1Height,wave1Dir,wave1Perturbation,ii*volumeTime*wave1Speed);
    [Uwave2,Vwave2,Wwave2] = waveGenerator(modelX,modelY,modelZ,wave2Length,wave2Height,wave2Dir,wave2Perturbation,ii*volumeTime*wave2Speed);
    [ interpolantsWave1 ] = createInterpolants( modelX,modelY,modelZ,Uwave1+Uwave2,Vwave1+Vwave2,Wwave1+Wwave2);
    [ DVwave1 ] = calculateDV( radarStruct, interpolantsWave1);
    [ DVwave2 ] = calculateDV( radarStruct_2, interpolantsWave1);
    
    maxW = ceil(max((Wwave1(:)+Wwave2(:)).*10))./10; %set scale max by rounding up to nearest tenth
    maxUV= ceil(max((DVwave1(:)).*10))./10;
    
    domainWidth = distance(radarLoc(2),waveEdgeW,radarLoc(2),waveEdgeE,referenceEllipsoid('earth','km'));
    
    if ii>0
        
        tt = ii.*(volumeTime/60);%for time in title
        
        FH=figure('Renderer', 'painters','Position',[10 10 600 600]);
        TLH = tiledlayout(FH,4,2,'Padding','compact','TileSpacing','compact');
        
        AH1 = nexttile(1,[2 1]);
        surf(squeeze(Xdist),squeeze(Ydist),squeeze(DVwave1),'EdgeColor','none','FaceColor','texturemap');view(2);
        set(gca,'box','on','XGrid','off','YGrid','off','ZGrid','off','FontSize',14,'xtick',-150:50:150,'ytick',-150:50:150,'Layer','top');
        ylabel('Distance (km)','FontSize',14); xlabel('Distance (km)','FontSize',14);
        colormap(gca,flipud(Colormap_DV_RdBu(64)));
        caxis([-maxUV maxUV]);colorbar;
        title(['Time: +' sprintf('%u:%02u',floor(tt),floor(((tt-floor(tt)).*60))) ' minutes'],'FontSize',14);
        axis image;
        hold on
        plot3([0 50],[0 0],[10 10],'LineWidth',2,'Color',[0.5 0.5 0.5]);
        
        
        AH2 = nexttile(2,[2 1]);
        surf(squeeze(Xdist),squeeze(Ydist),squeeze(DVwave1 - DVwaveLAST),'EdgeColor','none','FaceColor','texturemap');view(2);axis image;
        set(gca,'box','on','XGrid','off','YGrid','off','ZGrid','off','FontSize',14,'xtick',-150:50:150,'ytick',-150:50:150,'Layer','top');
        ylabel('Distance (km)','FontSize',14); xlabel('Distance (km)','FontSize',14);
        colormap(gca,[0 0 0;1 1 1]);colorbar('Visible','off');
        caxis([-0.065 -0.06]);
        title(['Time: +' sprintf('%u:%02u',floor(tt),floor(((tt-floor(tt)).*60))) ' minutes'],'FontSize',14);
        axis image;
        hold on
        plot3([0 50],[0 0],[10 10],'LineWidth',2,'Color',[0.5 0.5 0.5]);
        
        %AH3 = nexttile(5);
        %%AH1 = axes('Position',[1/16 6.5/8 5/16 1/8]);
        %imagesc([-1.*domainWidth./2 domainWidth./2],[0 domainHeight],squeeze(Wwave1(round(end/2),:,:)+Wwave2(round(end/2),:,:))');
        %%imagesc(squeeze(Xdist),[0 domainHeight],squeeze(Wwave1(300,:,:))');
        %set(AH3,'YDir','normal','FontSize',14);
        %ylabel('Height (m)','FontSize',14);
        %xlabel('Distance From Radar (km)','FontSize',14);
        %ylim([0 6000]);xlim([0 50]);
        %colormap(AH3,circshift(fliplr(Colormap_DV_RdBu(64)),2,2));
        %caxis([-maxW maxW]);colorbar;
        %title('Idealized Wave Simulated Vertical Velocity (m/s)','FontSize',14,'Interpreter','none');
        
        %AH4 = nexttile(7);
        AH4 = nexttile(5, [1 2]);
        %AH1 = axes('Position',[1/16 6.5/8 5/16 1/8]);
        imagesc([-1.*domainWidth./2 domainWidth./2],[0 domainHeight],squeeze(Uwave1(round(end/2),:,:)+Uwave2(round(end/2),:,:))');
        %imagesc(squeeze(Xdist),[0 domainHeight],squeeze(Wwave1(300,:,:))');
        set(AH4,'YDir','normal','FontSize',14);
        ylabel('Height (m)','FontSize',14);
        xlabel('Distance From Radar (km)','FontSize',14);
        ylim([0 6000]);xlim([0 100]);%xlim([0 50]);
        colormap(AH4,flipud(Colormap_DV_RdBu(64)));
        caxis([-maxUV maxUV]);colorbar;
        title('Idealized Wave Simulated Zonal Velocity (m/s)','FontSize',14,'Interpreter','none');
        
        %AH5 = nexttile(6);
        AH5 = nexttile(7);
        surf(squeeze(Xdist2(:,6,:)),squeeze(radarStruct_2.hRadar(:,6,:)),squeeze(DVwave2(:,6,:)),'EdgeColor','none','FaceColor','texturemap');view(2);
        set(gca,'box','on','XGrid','off','YGrid','off','ZGrid','off','FontSize',14,'xtick',0:25:150,'Layer','top');
        ylabel('Height (m)','FontSize',14); xlabel('Distance (km)','FontSize',14);
        colormap(gca,flipud(Colormap_DV_RdBu(64)));
        caxis([-maxUV maxUV]);colorbar;
        title(['Time: +' sprintf('%u:%02u',floor(tt),floor(((tt-floor(tt)).*60))) ' minutes'],'FontSize',14);
        hold on
        plot3([0 radarRange/1000],[wave1Height wave1Height],[10 10],'LineWidth',2,'Color',[0.5 0.5 0.5],'LineStyle', '--')
        plot3([0 radarRange/1000],[wave2Height wave2Height],[10 10],'LineWidth',2,'Color',[0.5 0.5 0.5],'LineStyle', '--')
        
        AH6 = nexttile(8);
        surf(squeeze(Xdist2(:,6,:)),squeeze(radarStruct_2.hRadar(:,6,:)),squeeze(DVwave2(:,6,:) - DVwaveLAST2(:,6,:)),'EdgeColor','none','FaceColor','texturemap');view(2);
        set(gca,'box','on','XGrid','off','YGrid','off','ZGrid','off','FontSize',14,'xtick',0:25:150,'Layer','top');
        ylabel('Height (m)','FontSize',14); xlabel('Distance (km)','FontSize',14);
        colormap(gca,[0 0 0;1 1 1]);colorbar('Visible','off');
        caxis([-0.065 -0.06]);
        title(['Time: +' sprintf('%u:%02u',floor(tt),floor(((tt-floor(tt)).*60))) ' minutes'],'FontSize',14);
        hold on
        plot3([0 radarRange/1000],[wave1Height wave1Height],[10 10],'LineWidth',2,'Color',[0.5 0.5 0.5],'LineStyle', '--')
        plot3([0 radarRange/1000],[wave2Height wave2Height],[10 10],'LineWidth',2,'Color',[0.5 0.5 0.5],'LineStyle', '--')
        
        if saveFigures
            FH.PaperUnits = 'inches';
            FH.PaperPosition = [0 0 12 12];
            print([savePath 'Demo_2_Wave_noBG_frame' num2str(ii,'%02u') '.png'],'-dpng','-r300');
            if ii == 1 && plotVector
                print([savePath 'Demo_2_Wave_noBG_frame' num2str(ii,'%02u') '.eps'],'-depsc');
                %exportgraphics(FH,[savePath 'Demo_2_Wave_noBG_frame' num2str(ii,'%02u') '.eps'],'ContentType','vector');
                %exportgraphics(FH,[savePath 'Demo_2_Wave_noBG_frame' num2str(ii,'%02u') '.pdf'],'ContentType','vector');
            end
        end
        
        close
    end
    DVwaveLAST = DVwave1;
    DVwaveLAST2 = DVwave2;
    
end


% vidBG=uint8(ones(1080,1920,3)).*255;
% writerObj = VideoWriter([savePath 'Demo_2_Wave_noBG_animation'],'MPEG-4');
% writerObj.FrameRate = 6;
% open(writerObj);
% frameList=dir([savePath 'Demo_2_Wave_noBG_frame*.png']);
% for ii=1:length(frameList)
%     Img=imread([savePath frameList(ii).name]);
%     IMG_resize = imresize(Img,[1080 NaN]);
%     [~,ImgX,~]=size(IMG_resize);
%     Img_out = vidBG;
%     colStart = 960-round((ImgX/2));
%     Img_out(:,colStart:colStart+ImgX-1,:) = IMG_resize;
%     writeVideo(writerObj,Img_out);
% end
% close(writerObj);

