function [DICpara] = setDICParas_STAQ(DICpara, file_name,Img,maskLeft,maskRight,LoadImgMethod)
%% ==============================================
% Decide DIC subset parameters
% numImages = size(file_name,2);



% Choose ROI
% -------------------
% fprintf('\n');
% disp('--- Define ROI corner points at the top-left and the bottom-right ---')
% ImgTemp = imread([file_name{2,1},'\',file_name{1,1}]); % Revised by Zach
% if size(ImgTemp,3) > 3
%     ImgTemp = ImgTemp(:,:,1:3); % Delete alpha channel
% end
% imshow(ImgTemp);
% title('Click top-left and the bottom-right corner points','fontweight','normal','fontsize',16);
% 
% gridx = zeros(1,2); gridy = zeros(1,2);
% [gridx(1), gridy(1)] = ginput(1);
% fprintf('Coordinates of top-left corner point are (%4.3f,%4.3f)\n',gridx(1), gridy(1))
% 
% [gridx(2), gridy(2)] = ginput(1);
% fprintf('Coordinates of bottom-right corner point are (%4.3f,%4.3f)\n',gridx(2), gridy(2))
% 
% gridxy.gridx = round(gridx); gridxy.gridy = round(gridy);
% ---------------------------
if DICpara.DICIncOrNot == 0 % acc
    [rows, cols] = find(maskLeft{1});
    gridx(1) = min(rows);
    gridy(1) = min(cols); 
    gridx(2) = max(rows);
    gridy(2) = max(cols); 
    gridxy.gridx = gridx; gridxy.gridy = gridy;

    % [rows_left, cols_left] = find(maskLeft{1});
    % [rows_right, cols_right] = find(maskRight{1});
    % gridx(1) = min([rows_left;rows_right]);
    % gridy(1) = min([cols_left;cols_right]); 
    % gridx(2) = max([rows_left;rows_right]);
    % gridy(2) = max([cols_left;cols_right]); 
    % gridxy.gridx = gridx; gridxy.gridy = gridy;

    % Check (Zach)
    % figure; imshow(maskLeft{1});
    % rectangle_matrix = zeros(size(maskLeft{1}));
    % if ~isempty(rows)
    %     rectangle_matrix( gridxy.gridx(1): gridxy.gridx(2), gridxy.gridy(1):gridxy.gridy(2)) = 1;
    % end
    %     figure; imshow(rectangle_matrix)
else
    gridxy.gridx = [1,size(maskLeft{1},1)]; gridxy.gridy = [1,size(maskLeft{1},2)];
end



% Choose subset size
fprintf('\n');
fprintf('--- What is the subset size? --- \n');
fprintf('Each subset has an area of [-winsize/2:winsize/2, -winsize/2:winsize/2] \n');
prompt = 'Input an even number: ';
winsize = input(prompt);

% Choose subset size
fprintf('--- What is the subset step? --- \n');
prompt = 'Input an integer to be a power of 2 (E.g., 16): ';
winstepsize = input(prompt);


% ==============================================
% Subproblem 2 solver: finite difference or finite element
Subpb2FDOrFEM = 0; % By default initialize parameters
Subpb2FDOrFEM = 1; % funParaInput('Subpb2FDOrFEM'); % Subproblem 2 using finite difference or fem?

% ==============================================
% Parallel cluster #
ClusterNo = funParaInput('ClusterNo'); % Assign parpool cluster No


% ==============================================
winsizeMin = funParaInput('winsizeMin'); % Assign the finest element size in the quadtree mesh
DICpara.winsizeMin = winsizeMin;


% ==============================================
% Deal with image sequence
% if numImages > 2
%     % ==============================================
%     % DIC initial guess
%     % NewFFTSearch = funParaInput('NewFFTSearch'); % Use last frame as init guess or not
% 
%     % ==============================================
%     % IncrementalOrNot  #
%     IncrementalOrNot = funParaInput('IncrementalOrNot');
%     % check if the requirements are met
%     if IncrementalOrNot ~= 0 % inc. mode
%         if (size(maskLeft,1) >= size(file_name,2))  &&  (size(MaskRight,1) >= size(file_name,2))
%         else
%             error('No enough mask files!');
%         end
%     else  % acc. mode: only need the first left and right mask
%         if (size(maskLeft,1) < 1)  &&  (size(MaskRight,1) < 1)
%             error('No enough mask files!');
%         end
%     end
% 
% 
%     try
%         switch IncrementalOrNot
%             case 0 % acc. default
%                 ImgSeqIncUnit = numImages+1;
%                 ImgSeqIncROIUpdateOrNot = 1;
%                 NewFFTSearch = 0;
%             case 1 % inc.
%                 ImgSeqIncUnit = 1;
%                 ImgSeqIncROIUpdateOrNot = 0; 
%                 NewFFTSearch = 1;
%             otherwise
%                 ImgSeqIncUnit = numImages+1;
%                 ImgSeqIncROIUpdateOrNot = 1;
%         end
% 
%     catch
%         ImgSeqIncUnit = numImages+1;
%         ImgSeqIncROIUpdateOrNot = 1;
%     end
% 
% 
% 
%     % ================================
% else % Only two frames
% 
%     ImgSeqIncUnit = numImages+1;
%     ImgSeqIncROIUpdateOrNot = 1;
%     IncrementalOrNot = 0; % Revised by Zach
% end

DICpara.winsize = winsize;
DICpara.winstepsize = winstepsize;
DICpara.gridxyROIRange = gridxy;
DICpara.LoadImgMethod = LoadImgMethod;
% DICpara.DICIncOrNot = IncrementalOrNot; % Revised by Zach
% DICpara.ImgSeqIncUnit = ImgSeqIncUnit;
% DICpara.ImgSeqIncROIUpdateOrNot = ImgSeqIncROIUpdateOrNot;
DICpara.Subpb2FDOrFEM = Subpb2FDOrFEM;
DICpara.ClusterNo = ClusterNo;
DICpara.ImgSize = size(Img{1});
DICpara.transformDisp = [];

end

