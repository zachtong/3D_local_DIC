function [U,F,HtempPar,LocalTime,ConvItPerEle,LocalICGNBadPtNum] = LocalICGN(U0,coordinatesFEM,...
    Df,ImgRef,ImgDef,DICpara,ICGNmethod,tol,shapeFunctionOrder)
%FUNCTION [U,F,HtempPar,LocalTime,ConvItPerEle,LocalICGNBadPtNum] = LocalICGN(U0,coordinatesFEM,...
%                                       Df,imgfNormalizedbc,imggNormalizedbc,DICpara,ICGNmethod,tol)
% The Local ICGN subset solver (part I): to assign a sequential computing or a parallel computing
% (see part II: ./func/funICGN.m)
% ----------------------------------------------
%   INPUT: U0                   Initial guess of the displacement fields
%          coordinatesFEM       FE mesh coordinates
%          Df                   Image grayscale value gradients
%          ImgRef               Reference image
%          ImgDef               Deformed image
%          DICpara              DIC parameters: subset size, subset spacing, ...
%          ICGNmethod           ICGN iteration scheme: 'GaussNewton' -or- 'LevenbergMarquardt'
%          tol                  ICGN iteration stopping threshold
%
%   OUTPUT: U                   Disp vector: [Ux_node1, Uy_node1, ... , Ux_nodeN, Uy_nodeN]';
%           F                   Deformation gradient tensor
%                               F = [F11_node1, F21_node1, F12_node1, F22_node1, ... , F11_nodeN, F21_nodeN, F12_nodeN, F22_nodeN]';
%           HtempPar            Hessian matrix for each local subset
%           LocalTime           Computation time
%           ConvItPerEle        ICGN iteration step for convergence
%           LocalICGNBadPtNum   Number of subsets whose ICGN iterations don't converge
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
% Last time updated: 02/2020.
% ==============================================


%% Initialization
warning('off');
winsize = DICpara.winsize;
winstepsize = DICpara.winstepsize;
ClusterNo = DICpara.ClusterNo;

temp = zeros(size(coordinatesFEM,1),1); UtempPar = temp; VtempPar = temp; % UtempPar = UPar{1}; VtempPar = UPar{2};
F11tempPar = temp; F21tempPar = temp; F12tempPar = temp; F22tempPar = temp;% F11tempPar = FPar{1}; F21tempPar = FPar{2}; F12tempPar = FPar{3}; F22tempPar = FPar{4};
if shapeFunctionOrder == 1
    HtempPar = zeros(size(coordinatesFEM,1),21);
elseif shapeFunctionOrder == 2
    HtempPar = zeros(size(coordinatesFEM,1),78);
end
ConvItPerEle = zeros(size(coordinatesFEM,1),1);

% -------- How to change parallel pools ---------
% myCluster = parcluster('local');
% myCluster.NumWorkers = 4;  % 'Modified' property now TRUE
% saveProfile(myCluster);    % 'local' profile now updated,
%                            % 'Modified' property now FALSE
% -------------- Or we can do this --------------
% Go to the Parallel menu, then select Manage Cluster Profiles.
% Select the "local" profile, and change NumWorkers to 4.
% -----------------------------------------------

%% ClusterNo == 0 or 1: Sequential computing
if (ClusterNo == 0) || (ClusterNo == 1)


    h = waitbar(0,'Please wait for Subproblem 1 IC-GN iterations!'); tic;

    for tempj = 1:size(coordinatesFEM,1)  % tempj is the element index
        try
            x0temp = coordinatesFEM(tempj,1); y0temp = coordinatesFEM(tempj,2);

            if shapeFunctionOrder == 1 % 1st-order shape function
                    [Utemp, Ftemp, ConvItPerEle(tempj), HtempPar(tempj,:)] = funICGN(...
                        U0(2*tempj-1:2*tempj),x0temp,y0temp,Df,ImgRef,ImgDef,winsize,tol,ICGNmethod);

            elseif shapeFunctionOrder ==2 % 2nd-order shape function
                [Utemp, Ftemp, ConvItPerEle(tempj), HtempPar(tempj,:)] = funICGN2(...
                    U0(2*tempj-1:2*tempj),x0temp,y0temp,Df,ImgRef,ImgDef,winsize,tol,ICGNmethod);

            end
            % disp(['ele ',num2str(tempj),' converge step is ',num2str(ConvItPerEle(tempj)),' (>0-converged; 0-unconverged)']);

            % ------ Store solved deformation gradients ------
            if shapeFunctionOrder == 1 % 2nd-order doesn't need F
                UtempPar(tempj) = Utemp(1); VtempPar(tempj) = Utemp(2);
                F11tempPar(tempj) = Ftemp(1); F21tempPar(tempj) = Ftemp(2); F12tempPar(tempj) = Ftemp(3); F22tempPar(tempj) = Ftemp(4);
                waitbar(tempj/(size(coordinatesFEM,1)));
            elseif shapeFunctionOrder == 2 % TBD
                UtempPar(tempj) = Utemp(1); VtempPar(tempj) = Utemp(2);
                F11tempPar(tempj) = Ftemp(1); F21tempPar(tempj) = Ftemp(2); F12tempPar(tempj) = Ftemp(3); F22tempPar(tempj) = Ftemp(4);
                % TBD
            end
        catch
            ConvItPerEle(tempj) = -1;
            UtempPar(tempj) = nan; VtempPar(tempj) = nan;
            F11tempPar(tempj) = nan; F21tempPar(tempj) = nan; F12tempPar(tempj) = nan; F22tempPar(tempj) = nan;
            waitbar(tempj/(size(coordinatesFEM,1)));
        end
    end
    close(h); LocalTime = toc;

    %% ClusterNo > 1: parallel computing
else

    % Start parallel computing
    % ****** This step needs to be careful: may be out of memory ******
    disp('--- Set up Parallel pool ---'); tic;
    hbar = parfor_progressbar(size(coordinatesFEM,1),'Please wait for Subproblem 1 IC-GN iterations!');
    
    tempU0 = [U0(1:2:end),U0(2:2:end)];

    for tempj = 1:size(coordinatesFEM,1)  % tempj is the element index
        try
            % Build Interp. coef. matrix.
            % For parallel running, all pixels has their own CoefMatrix.
            % BicubicBspline = BicubiclBsplineInit(ImgDef);

            tempcoordinatesFEM = coordinatesFEM(tempj,:);
            x0temp = tempcoordinatesFEM(1); y0temp = tempcoordinatesFEM(2);
            if shapeFunctionOrder == 1
                % Bicubic Bspline
                [Utemp, Ftemp, ConvItPerEle(tempj), HtempPar(tempj,:)] = funICGN(...
                    tempU0(tempj,:),x0temp,y0temp,Df,ImgRef,ImgDef,winsize,tol,ICGNmethod);

           
            elseif shapeFunctionOrder == 2
                % Since 2nd-order shape function is just used for stereo
                % matching for now. We don't need Ftemp.
                [Utemp, Ftemp, ConvItPerEle(tempj), HtempPar(tempj,:)] = funICGN2(U0(2*tempj-1:2*tempj), ...
                    x0temp,y0temp,Df,ImgRef,ImgDef,winsize,tol,ICGNmethod,BicubicBsplineConst.Value);

                % Old code by Zach
                % Create interpolation mapping
                % InterpCoef=griddedInterpolant({1:1:size(ImgDef',1),1:1:size(ImgDef',2)}, ImgDef', 'spline');
                % [Utemp, Ftemp, ConvItPerEle(tempj), HtempPar(tempj,:)] = funICGN2_backup_spline(U0(2*tempj-1:2*tempj), ...
                % x0temp,y0temp,Df,ImgRef,ImgDef,winsize,tol,ICGNmethod);
            end
            % disp(['ele ',num2str(tempj),' converge step is ',num2str(ConvItPerEle(tempj)),' (>0-converged; 0-unconverged)']);
            % ------ Store solved deformation gradients ------
            UtempPar(tempj) = Utemp(1); VtempPar(tempj) = Utemp(2);
            if shapeFunctionOrder == 1 % 2nd-order doesn't need F
                F11tempPar(tempj) = Ftemp(1); F21tempPar(tempj) = Ftemp(2); F12tempPar(tempj) = Ftemp(3); F22tempPar(tempj) = Ftemp(4);
                hbar.iterate(1);
            elseif shapeFunctionOrder == 2 % TBD
                % F11tempPar(tempj) = Ftemp(1); F21tempPar(tempj) = Ftemp(2); F12tempPar(tempj) = Ftemp(3); F22tempPar(tempj) = Ftemp(4);
                % hbar.iterate(1);
            end
        catch ME
            fprintf('Error ID: %s\n', ME.identifier);
            fprintf('Error Message: %s\n', ME.message);
            ConvItPerEle(tempj) = -1;
            UtempPar(tempj) = nan; VtempPar(tempj) = nan;
            F11tempPar(tempj) = nan; F21tempPar(tempj) = nan; F12tempPar(tempj) = nan; F22tempPar(tempj) = nan;
            hbar.iterate(1);
        end
    end

    close(hbar);
    LocalTime = toc;

end

U = U0; U(1:2:end) = UtempPar; U(2:2:end) = VtempPar;

if shapeFunctionOrder == 1
    F = zeros(4*size(coordinatesFEM,1),1); F(1:4:end) = F11tempPar; F(2:4:end) = F21tempPar; F(3:4:end) = F12tempPar; F(4:4:end) = F22tempPar;
    % ------ Clear bad points for Local DIC ------
    % find bad points after Local Subset ICGN
    [row1,~] = find(ConvItPerEle(:)<0);
    [row2,~] = find(ConvItPerEle(:)>99);
    [row3,~] = find(ConvItPerEle(:)==102);
    LocalICGNBadPt = unique(union(row1,row2)); LocalICGNBadPtNum = length(LocalICGNBadPt)-length(row3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Though some subsets are converged, but their accuracy is worse than most
    % other subsets. This step is to remove those subsets with abnormal convergence steps
    LocalICGNGoodPt = setdiff([1:1:size(coordinatesFEM,1)],LocalICGNBadPt);
    ConvItPerEleMean = mean(ConvItPerEle(LocalICGNGoodPt));
    ConvItPerEleStd = std(ConvItPerEle(LocalICGNGoodPt));
    [row4,~] = find(ConvItPerEle(:) > max([ConvItPerEleMean+0.15*ConvItPerEleStd ])); % Here "0.15" is an empirical value
    LocalICGNBadPt = unique(union(LocalICGNBadPt,row4));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Print results info on the MATLAB command window
    disp(['Local ICGN bad subsets %: ', num2str(LocalICGNBadPtNum),'/',num2str(size(coordinatesFEM,1)-length(row3)), ...
        '=',num2str(100*(LocalICGNBadPtNum)/(size(coordinatesFEM,1)-length(row3))),'%']);
    U(2*LocalICGNBadPt-1) = NaN; U(2*LocalICGNBadPt) = NaN;
    F(4*LocalICGNBadPt-3) = NaN; F(4*LocalICGNBadPt-2) = NaN; F(4*LocalICGNBadPt-1) = NaN; F(4*LocalICGNBadPt) = NaN;
    % figure, scatter(coordinatesFEM(:,1),coordinatesFEM(:,2),[],U(1:2:end));
    % figure, scatter(coordinatesFEM(:,1),coordinatesFEM(:,2),[],F(4:4:end));
elseif shapeFunctionOrder == 2
    % TBD
    % Since 2nd-order shape function is just used for stereo matching for
    % now. We don't need Ftemp.
end

% ------ inpaint nans using gridfit ------
nanindex = find(isnan(U(1:2:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),U(2*notnanindex-1),'natural');
U1 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),U(2*notnanindex),'natural');
U2 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));

U = [U1(:),U2(:)]'; U = U(:);


if shapeFunctionOrder == 1
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),F(4*notnanindex-3),'natural');
    F11 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),F(4*notnanindex-2),'natural');
    F21 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),F(4*notnanindex-1),'natural');
    F12 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),F(4*notnanindex-0),'natural');
    F22 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));

    F = [F11(:),F21(:),F12(:),F22(:)]'; F = F(:);
elseif shapeFunctionOrder == 2 % TBD
    F = []; LocalICGNBadPtNum = [];
end

end






