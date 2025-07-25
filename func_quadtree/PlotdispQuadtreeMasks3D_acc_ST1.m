function PlotdispQuadtreeMasks3D_acc_ST1(U_3D,U_2D_L,FirstFEM,FirstImg,CurrentImg,DICpara,voidIndex)
%PLOTDISPQUADTREE: to plot DIC solved displacement components  
%   PlotdispQuadtree(U,coordinatesFEMWorld,elementsFEM,CurrentImg,DICpara)
% ----------------------------------------------
%
%   INPUT: U                    Displacement vector: 
%                               U = [Ux_node1, Uy_node1, Ux_node2, Uy_node2, ... , Ux_nodeN, Uy_nodeN]';
%          coordinatesFEM       FE mesh coordinates
%          elementsFEM          FE mesh elements
%          CurrentImg           Current deformed image
%          DICpara              DIC paramters
%
%   OUTPUT: Plots of x-displacement field and y-displacement field.
%
% ----------------------------------------------
% Reference
% [1] RegularizeNd. Matlab File Exchange open source. 
% https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
% [2] Gridfit. Matlab File Exchange open source. 
% https://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last date modified: 2020.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');
% run('./plotFiles/Black_rainbow.m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% convert pixel unit to the physical world unit %%%%%
try um2px = DICpara.um2px; 
catch um2px = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OrigDICImgTransparency = DICpara.OrigDICImgTransparency; % Original raw DIC image transparency

Image2PlotResults = DICpara.Image2PlotResults; % Choose image to plot over (first only, second and next images)

%%%% Important !!!!!%%%%%%%%%%%
% Currently we don't support plot on the deformed images (Zach 2024.3.5)
% 原因：目前我算的位移全都是基于初始frame的，如果要在其他的frame展示，由于每一个frame的MESH不一样，
% SHOW函数要求我们把位移也全部先插值到其他的frame上，会增加很多运算成本
% PS. 应该有不用插值的方法，需要在temporal matching的时候，把2D_inc_disp结果一起输出出来，然后在
% 3D重构的时候用2d_inc_disp/coor 得到3d_inc_coor，然后再计算3d_inc_disp（仍然需要转换出3d_acc_disp）。
% 此时的3d_inc_disp就可以用于这个函数和plot_strain函数的show函数。



%% Plot on deformed images or Not
if Image2PlotResults == 1
    coordinatesFEMWorldDef = [FirstFEM.coordinatesFEM(:,1) + U_2D_L(:,1), FirstFEM.coordinatesFEM(:,2) + U_2D_L(:,2)];
    elementsFEM = FirstFEM.elementsFEM;
    Img = CurrentImg;
    disp_u = U_3D{1}; disp_v = U_3D{2}; disp_w = U_3D{3};
else
    coordinatesFEMWorldDef = [FirstFEM.coordinatesFEM(:,1), FirstFEM.coordinatesFEM(:,2)];
    elementsFEM = FirstFEM.elementsFEM;
    Img = FirstImg;
    disp_u = U_3D{1}; disp_v = U_3D{2}; disp_w = U_3D{3};
end

%%%%%%%%%%% JY!!!Mask START %%%%%%%%%%%%%%%
% if Image2PlotResults == 1
%     if ~isempty(CurrentImgMask)
%         for tempi = 1:size(coordinatesFEMWorldDef,1)
%             try
%                 if CurrentImgMask( round(coordinatesFEMWorldDef(tempi,1)/um2px), ...
%                                     (size(CurrentImgMask,2)+1-round(coordinatesFEMWorldDef(tempi,2)/um2px)) ) == 0 
%                     coordinatesFEMWorldDef(tempi,:) = [nan,nan];
%                 end
%             catch
%                 coordinatesFEMWorldDef(tempi,:) = [nan,nan];
%             end
% 
%         end
%     else
%         CurrentImgMask = imread(CurrentImg)';
%         for tempi = 1:size(coordinatesFEMWorldDef,1)
%             try
%                 if CurrentImgMask( round(coordinatesFEMWorldDef(tempi,1)/um2px), ...
%                         (size(CurrentImgMask,2)+1-round(coordinatesFEMWorldDef(tempi,2)/um2px)) ) < 0
%                     coordinatesFEMWorldDef(tempi,:) = [nan,nan];
%                 end
%             catch
%                 coordinatesFEMWorldDef(tempi,:) = [nan,nan];
%             end
%         end
%     end
% end
%%%%%%%%%%% JY!!!Mask END %%%%%%%%%%%%%%%



% Zach: Warpage = w - w_mean
%disp_w = disp_w - mean(disp_w) ;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 1) dispx u ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes;  
try h1=imshow( flipud(imread(Img)) ,'InitialMagnification','fit');
catch h1=surf(  flipud( imread(Img) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([], elementsFEM(:,1:4), coordinatesFEMWorldDef/um2px, disp_u, 'NoEdgeColor');
set(gca,'fontSize',18); set(gca,'ydir','reverse');view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency);  colormap("turbo"); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(cMap); caxis([-40,40]); % colormap(jet);  
% caxis([-35,35]); % caxis([-0.025,0.025]); 
% caxis([-10,10]); D shaped
%caxis([-0.08,0.03]);  % Bulge
%caxis([-2.5 0]);  % Pig heart

 %caxis([-1 0.5]);
% colormap(black_rainbow);  
%  colormap(jet); caxis([-20 20]);
% ax1.XTick = [100,200,300]; % Unit: px
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; % ax1.TickLabelInterpreter = 'latex'; 

%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick))');
cb2 = colorbar('Position',[.17+0.685+0.012 .11+.128 .03 .557 ]); % cb2.TickLabelInterpreter = 'latex';
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 2) dispy v ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes;  
try h1=imshow( flipud(imread(Img)) ,'InitialMagnification','fit');
catch h1=surf(  flipud( imread(Img) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,disp_v,'NoEdgeColor');
set(gca,'fontSize',18);set(gca,'ydir','reverse'); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency);  colormap("turbo"); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(cMap); caxis([-40 ,40]); % colormap(jet);  
% caxis([-0.7,0.7]);
% caxis([0,10]); D shaped
%caxis([-0.1 0.05]) % bulge
%caxis([-2 2]) % Pig heart

% caxis([-0.7,0.7]); % caxis([-0.025,0.025]); 
% colormap(black_rainbow);    caxis([-0.5,0]);
%   colormap(jet); caxis([-20 20]);
% ax1.XTick = [100,200,300]; % Unit: px
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; % ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
cb2 = colorbar('Position',[.17+0.685+0.012 .11+.128 .03 .557 ]); % cb2.TickLabelInterpreter = 'latex';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 3) dispz w ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes;  
try h1=imshow( flipud(imread(Img)) ,'InitialMagnification','fit');
catch h1=surf(  flipud( imread(Img) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,disp_w,'NoEdgeColor');
set(gca,'fontSize',18);set(gca,'ydir','reverse'); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency);  colormap("turbo"); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(cMap); caxis([-40 ,40]); % colormap(jet);  
% caxis([-0.7,0.7]);
% caxis([-0.3,0.3]); D shaped
%caxis([0 1.2]) %bulge
%caxis([-1 3]) % pig heart
% caxis([-0.7,0.7]); % caxis([-0.025,0.025]); 
% colormap(black_rainbow);    caxis([-0.5,0]);
%   colormap(jet); caxis([-20 20]);
% ax1.XTick = [100,200,300]; % Unit: px
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; % ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
cb2 = colorbar('Position',[.17+0.685+0.012 .11+.128 .03 .557 ]); % cb2.TickLabelInterpreter = 'latex';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 4) disp magnitude ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig1=figure; ax1=axes;  
% try h1=imshow( flipud(imread(Img)) ,'InitialMagnification','fit');
% catch h1=surf(  flipud( imread(Img) ),'EdgeColor','none','LineStyle','none');
% end
% 
% axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
% hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,sqrt(disp_u.^2+disp_v.^2+disp_w.^2),'NoEdgeColor');
% set(gca,'fontSize',18);set(gca,'ydir','reverse'); view(2); box on; axis equal;  axis tight;   
% alpha(h2,OrigDICImgTransparency);  colormap(turbo); caxis auto;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% TODO: manually modify colormap and caxis %%%%%%
% % colormap(cMap); caxis([-40 ,40]); % colormap(jet);  
% % caxis([-0.7,0.7]);
% % caxis([-0.7,0.7]); % caxis([-0.025,0.025]); 
% % colormap(black_rainbow);    caxis([-0.5,0]);
% %   colormap(jet); caxis([-20 20]);
% % ax1.XTick = [100,200,300]; % Unit: px
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on'; % ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
% % cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
% cb2 = colorbar('Position',[.17+0.685+0.012 .11+.128 .03 .557 ]); % cb2.TickLabelInterpreter = 'latex';



