%%% INPUTS %%%

%I/O Settings
%savePath = 'D:\Dropbox\matlab\RadarSim\WaveDemo\WaveDepth\';

%Radar Parameters
Tilts=0.5:0.25:35;
[ radarStruct ] = createRadarStructure(-72.859, 40.861, 15, Tilts, 80000, 200, 35.9);%KOKX
[Xdist,Xaz] = distance(40.861,-72.859,40.861,radarStruct.lonRadar,referenceEllipsoid('earth','km'));Xdist(Xaz>180)=Xdist(Xaz>180).*-1;
[Ydist,Yaz] = distance(40.861,-72.859,radarStruct.latRadar,-72.859,referenceEllipsoid('earth','km'));Ydist(Yaz>90)=Ydist(Yaz>90).*-1;

%Wave Parameters
WL1=6000;  %wavelength (m)
WA = 3000;  %wave altitude (m)
WT = 1000;  %wave thickness (m)
Dir1 = 90;  %wave direction (deg cw from north)
WSpeed1=15;  %wave propegation speed (m/s)
%WSpeed2=8;  %wave propegation speed  (m/s)
WPurt1=3;  %wave horizontal velocity perturbation (m/s)
volumeTime=4*60;  %time between radar volumes (s)
bgRange = [-2 3];  %background wind parameters (m/s)

%%% END INPUTS %%%

[modelX,modelY,modelZ] = meshgrid(linspace(-73,-71.8,600),linspace(39.2,42.5,100),linspace(0,5000,225));


%make the waves
[Uwave1,Vwave1,Wwave1] = KHGenerator(modelX,modelY,modelZ,WL1,WA,WT,Dir1,WPurt1,bgRange,0);%,ii.*(60*10*(4+1/6)));
[ interpolantsWave1 ] = createInterpolants( modelX,modelY,modelZ,Uwave1,Vwave1,Wwave1);
[ DVwave1 ] = calculateDV( radarStruct, interpolantsWave1);

[Uwave1,Vwave1,Wwave1] = KHGenerator(modelX,modelY,modelZ,WL1,WA,WT,Dir1,WPurt1,bgRange,volumeTime*WSpeed1);%,ii.*(60*10*(4+1/6)));
[ interpolantsWave1 ] = createInterpolants( modelX,modelY,modelZ,Uwave1,Vwave1,Wwave1);
[ DVwave2 ] = calculateDV( radarStruct, interpolantsWave1);


%plotting
FH=figure;
subplot(3,1,1);
surf(abs(squeeze(Xdist(:,7,:))),squeeze(radarStruct.hRadar(:,7,:)),zeros(size(squeeze(radarStruct.hRadar(:,7,:)))),squeeze(DVwave1(:,7,:)),'EdgeColor','none','FaceColor','texturemap');
set(gca,'box','on','XGrid','off','YGrid','off','ZGrid','off','FontSize',16,'Layer','top','DataAspectRatio', [1 1000 1]);
view(2)
colormap(gca,flipud(Colormap_DV_RdBu(64)));
caxis([-7 7]);colorbar;
hold on
ylim([0 5000]);xlim([0 40]);
xlabel('Distance From Radar (km)','FontSize',16);
ylabel('Height (m)','FontSize',16);
title('KH Radial Velocity Cross Section','FontSize',20,'Interpreter','none');

subplot(3,1,2);
surf(abs(squeeze(Xdist(:,7,:))),squeeze(radarStruct.hRadar(:,7,:)),zeros(size(squeeze(radarStruct.hRadar(:,7,:)))),squeeze(DVwave2(:,7,:)-DVwave1(:,7,:)),'EdgeColor','none','FaceColor','texturemap');
set(gca,'box','on','XGrid','off','YGrid','off','ZGrid','off','FontSize',16,'Layer','top','DataAspectRatio', [1 1000 1]);
view(2)
colormap(gca,flipud(Colormap_DV_RdBu(64)));
caxis([-7 7]);colorbar;
hold on
ylim([0 5000]);xlim([0 40]);
xlabel('Distance From Radar (km)','FontSize',16);
ylabel('Height (m)','FontSize',16);
title('KH Wave Radial Velocity Difference Cross Section','FontSize',20,'Interpreter','none');

subplot(3,1,3);
surf(abs(squeeze(Xdist(:,7,:))),squeeze(radarStruct.hRadar(:,7,:)),zeros(size(squeeze(radarStruct.hRadar(:,7,:)))),squeeze(DVwave2(:,7,:)-DVwave1(:,7,:)),'EdgeColor','none','FaceColor','texturemap');
set(gca,'box','on','XGrid','off','YGrid','off','ZGrid','off','FontSize',16,'Layer','top','DataAspectRatio', [1 1000 1]);
view(2)
colormap(gca,[0 0 0; 1 1 1]);colorbar('Visible','off');
caxis([-0.065 -0.060]);
hold on
ylim([0 5000]);xlim([0 40]);
xlabel('Distance From Radar (km)','FontSize',16);
ylabel('Height (m)','FontSize',16);
title('KH Wave Detection Cross Section','FontSize',20,'Interpreter','none');

%FH.PaperUnits = 'inches';
%FH.PaperPosition = [0 0 16 9];
%print([savePath 'KH_Demo.png'],'-dpng','-r120');
