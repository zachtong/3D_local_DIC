function [DICpara,x0,y0,u,v,cc]= IntegerSearchQuadtree(ImgRef,ImgDef,file_name,DICpara,ImgSeqNum)
%FUNCTION [DICpara,x0,y0,u,v,cc]= IntegerSearch(fNormalized,gNormalized,file_name,DICpara)
% Objective: To compute an inititial guess of the unknown displacement
% field by maximizing the FFT-based cross correlation
% ----------------------------------------------
%   INPUT: ImgRef       Reference image
%          ImgDef       Deformed image
%          file_name    Loaded DIC raw images file name
%          DICpara      Current DIC parameters
%
%   OUTPUT: DICpara     Updated DIC parameters
%           x0,y0       DIC subset x- and y- positions
%           u,v         x- and y- displacements
%           cc          Cross correlation information
%
% ----------------------------------------------
% Author: Jin Yang.
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 02/2020.
% ==============================================


%% Initialization
gridxROIRange = DICpara.gridxyROIRange.gridx;
gridyROIRange = DICpara.gridxyROIRange.gridy;
winsize = DICpara.winsize;
winstepsize = DICpara.winstepsize;
showImgOrNot = DICpara.showImgOrNot;

% Zach:
% switch camera0OrNot
%     case 'camera0'
%     ImgRefMask = DICpara.ImgRefMask_left;
%     case 'notCamera0'
%     ImgRefMask = DICpara.ImgRefMask_right;
% end
% figure;imshow(ImgRefMask);

try
    InitFFTSearchMethod = DICpara.InitFFTSearchMethod;
    % disp(num2str(InitFFTSearchMethod));
catch
    InitFFTSearchMethod = funParaInput('InitFFTSearchMethod');
end


%% To compute the inititial guess from maximizing the FFT-based cross correlation
if (InitFFTSearchMethod == 1) || (InitFFTSearchMethod == 2)
    InitialGuessSatisfied = 1;
    while InitialGuessSatisfied == 1

        fprintf('--- The size of initial guess search zone (pixels)? ---  \n')
        fprintf('User should start to try a small integer value, and gradually increase the value of \n');
        fprintf('the search zone size until it is larger than the magnitudes of |disp u| and |disp v|. \n');
        fprintf('User could also input [size_x, size_y] to search in a rectangular zone. \n');
        prompt = 'Input here: ';

        if ImgSeqNum == 2
            tempSizeOfSearchRegion = input(prompt);
            fprintf('Do you want to use fixed search distance? 0-yes 1-no \n');
            DICpara.fixSearchDistanceOrNot = input(prompt);
            %DICpara.fixSearchDistanceOrNot = 0;
            DICpara.NewFFTSearchDistance = tempSizeOfSearchRegion;

        elseif ImgSeqNum > 2 && DICpara.fixSearchDistanceOrNot == 1
            fprintf('Search distance in a rectangular zone. \n');
            tempSizeOfSearchRegion = input(prompt);
        elseif  ImgSeqNum > 2 && DICpara.fixSearchDistanceOrNot == 0
            tempSizeOfSearchRegion = DICpara.NewFFTSearchDistance;
        elseif ImgSeqNum == -1 % Stereomatching integer search
            tempSizeOfSearchRegion = input(prompt);
        end

        if length(tempSizeOfSearchRegion) == 1, tempSizeOfSearchRegion = tempSizeOfSearchRegion*[1,1]; end

        if (InitFFTSearchMethod == 1) % whole field for initial guess,
            [x0,y0,u,v,cc] = funIntegerSearch(ImgRef,ImgDef,tempSizeOfSearchRegion,gridxROIRange,gridyROIRange,winsize,winstepsize,0,winstepsize);

        else % (InitFFTSearchMethod == 2), several local seeds for initial guess

            % Open DIC image, and manually click several local seeds.
            figure; imshow( (imread(file_name{1})) ); % surf(fNormalized,'EdgeColor','none','LineStyle','none'); view(2);
            [row1, col1] = ginput; row = floor(col1); col = floor(row1);
            [x0,y0,u,v,cc] = funIntegerSearch(ImgRef,ImgDef,tempSizeOfSearchRegion,gridxROIRange,gridyROIRange,winsize,winstepsize,1,[row,col]);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Apply ImgRefMask to make u,v nans if there is a hole
        try
            x0y0Ind = sub2ind(DICpara.ImgSize, x0(:), y0(:));
            temp1 = double(DICpara.ImgRefMask(x0y0Ind));
            temp1(~logical(temp1))=nan;
            HolePtIndMat=reshape(temp1,size(x0));
            u = u.*HolePtIndMat; v = v.*HolePtIndMat;
        catch

        end
        % --------------------------------------
        if showImgOrNot == 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Have a look at the integer search results
            % --------------------------------------
            close all;
            figure; surf(u); colorbar;
            title('Displacement u','fontweight','normal')
            set(gca,'fontSize',18);
            title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
            axis tight; %axis equal; % set(gca,'XTick',[] );
            xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            set(gcf,'color','w');
            a = gca; a.TickLabelInterpreter = 'latex';
            b = colorbar; b.TickLabelInterpreter = 'latex';
            box on; colormap jet;
            pause(0.5);

            figure; surf(v); colorbar;
            title('Displacement v','fontweight','normal')
            set(gca,'fontSize',18);
            title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
            axis tight; %axis equal; % set(gca,'XTick',[] );
            xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            set(gcf,'color','w');
            a = gca; a.TickLabelInterpreter = 'latex';
            b = colorbar; b.TickLabelInterpreter = 'latex';
            box on; colormap jet;
            pause(0.5);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %        fprintf('--- Are you satisfied with initial guess with current search region? (0-yes; 1-no)? ---  \n')
            %        prompt = 'Input here: ';
            %        InitialGuessSatisfied = input(prompt);
            InitialGuessSatisfied = 0;
        else
            InitialGuessSatisfied = 0;

        end
    end

    % ======== Find some bad inital guess points ========
    cc.ccThreshold = 1.25; % bad cross-correlation threshold (mean - ccThreshold*stdev for q-factor distribution)
    qDICOrNot = 0;
    Thr0 = 100; [u,v,cc] = funRemoveOutliers(u,v,cc,qDICOrNot,Thr0);

    %%
elseif InitFFTSearchMethod == 0 % Multigrid search

    tempSizeOfSearchRegion = 0;
    [x0,y0,u,v,cc] = funIntegerSearchMg(ImgRef,ImgDef,gridxROIRange,gridyROIRange,winsize,winstepsize,winstepsize);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply ImgRefMask to make u,v nans if there is a hole
    try
        x0y0Ind = sub2ind(DICpara.ImgSize, x0(:), y0(:));
        temp1 = double(DICpara.ImgRefMask(x0y0Ind));
        temp1(~logical(temp1))=nan;
        HolePtIndMat=reshape(temp1,size(x0));
        u = u.*HolePtIndMat; v = v.*HolePtIndMat;
    catch
    end
    % --------------------------------------
    if showImgOrNot == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting initial guess
        % --------------------------------------
        close all;
        figure; surf(u); colorbar;
        title('Displacement u','fontweight','normal')
        set(gca,'fontSize',18);
        title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
        axis tight; %axis equal; % set(gca,'XTick',[] );
        xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
        set(gcf,'color','w');
        a = gca; a.TickLabelInterpreter = 'latex';
        b = colorbar; b.TickLabelInterpreter = 'latex';
        box on; colormap jet;
        pause(0.5);


        figure; surf(v); colorbar;
        title('Displacement v','fontweight','normal')
        set(gca,'fontSize',18);
        title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
        axis tight; %axis equal; % set(gca,'XTick',[] );
        xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
        set(gcf,'color','w');
        a = gca; a.TickLabelInterpreter = 'latex';
        b = colorbar; b.TickLabelInterpreter = 'latex';
        box on; colormap jet;
        pause(0.5);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
    end
elseif InitFFTSearchMethod == 3 % All initial guess are zero.
    x0 = DICpara.gridxyROIRange.gridx(1) : DICpara.winstepsize : DICpara.gridxyROIRange.gridx(2);
    x0 = repmat(x0, numel(DICpara.gridxyROIRange.gridy(1) : DICpara.winstepsize : DICpara.gridxyROIRange.gridy(2)), 1);

    y0 = (DICpara.gridxyROIRange.gridy(1):DICpara.winstepsize:DICpara.gridxyROIRange.gridy(2))';
    y0 = repmat(y0, 1, numel(DICpara.gridxyROIRange.gridx(1):DICpara.winstepsize:DICpara.gridxyROIRange.gridx(2)));

    u = zeros([size(x0,1),size(x0,2)]); v = zeros([size(x0,1),size(x0,2)]);
    cc.max = ones([size(x0,1),size(x0,2)]);
    tempSizeOfSearchRegion = [0,0];

end

% ======== Finally, update DICpara ========
DICpara.InitFFTSearchMethod = InitFFTSearchMethod;
DICpara.SizeOfFFTSearchRegion = tempSizeOfSearchRegion;


