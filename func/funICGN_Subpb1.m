function [U,H,stepwithinwhile,varargout] = funICGN_Subpb1(x0,y0,Df,ImgRef,ImgDef,winsize,...
    H,beta,mu,udual,vdual,UOld,FOld,tol,ICGNmethod,varargin)
%FUNCTION [U,H,stepwithinwhile] = funICGN_Subpb1(x0,y0,Df,ImgRef,ImgDef,winsize,...
%                                                H,beta,mu,udual,vdual,UOld,FOld,tol,ICGNmethod)
% The ALDIC Subproblem 1 ICGN subset solver (part II): ICGN iteration
% (see part I: ./func/Subpb1.m -or- ./func_quadtree/Subpb1Quadtree.m)
% ----------------------------------------------
%   INPUT: x0,y0                FE mesh nodal coordinates
%          Df                   Image grayscale value gradients
%          ImgRef               Reference image
%          ImgDef               Deformed image
%          winsize              DIC parameter subset size
%          H                    Stored Hessian matrix
%          mu,beta              ALDIC coefficients
%          udual,vdual          Dual variables
%          UOld                 Initial guess of the displacement fields
%          FOld                 Initial guess of the deformation gradients
%          tol                  ICGN iteration stopping threshold
%          ICGNmethod           ICGN iteration scheme: 'GaussNewton' -or- 'LevenbergMarquardt'
%
%   OUTPUT: U                   Disp vector: [Ux_node1, Uy_node1, ... , Ux_nodeN, Uy_nodeN]';
%           H                   Hessian matrix for each local subset
%           stepwithinwhile     ICGN iteration step for convergence
%
% ----------------------------------------------
% Author: Jin Yang.
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 2018.03, 2020.12.
% ==============================================


%% Initialization
warning('off');
DfCropWidth = Df.DfCropWidth;
imgSize = Df.imgSize;


%% ---------------------------
% Find local subset region
x  = [x0-winsize/2 ; x0+winsize/2 ; x0+winsize/2 ; x0-winsize/2];
y  = [y0-winsize/2 ; y0-winsize/2 ; y0+winsize/2 ; y0+winsize/2];

% ---------------------------
% Initialization: Get P0
P0 = [FOld(1) FOld(2) FOld(3) FOld(4) UOld(1) UOld(2)]'; P = P0;

% ---------------------------
% Find region for f
[XX,YY] = ndgrid([x(1):1:x(3)],[y(1):1:y(3)]);
%tempf = imgfNormalizedbc.eval(XX,YY);
%DfDx = imgfNormalizedbc.eval_Dx(XX,YY);
%DfDy = imgfNormalizedbc.eval_Dy(XX,YY);
tempf = ImgRef([x(1):1:x(3)],[y(1):1:y(3)]);
DfDx = Df.DfDx((x(1)-DfCropWidth):1:(x(3)-DfCropWidth), (y(1)-DfCropWidth):1:(y(3)-DfCropWidth));
DfDy = Df.DfDy((x(1)-DfCropWidth):1:(x(3)-DfCropWidth), (y(1)-DfCropWidth):1:(y(3)-DfCropWidth));


%% %%%%%%%% If there are >50% of the subset are painted with patterns %%%%%%%%%%%%
[DfDxImgMaskIndRow,~] = find(DfDx==0);
if length(DfDxImgMaskIndRow)<0.50*(winsize+1)^2

    if length(DfDxImgMaskIndRow)>0.1*(winsize+1)^2
        winsize = max(round(sqrt(2)*winsize));
        x = [x0-winsize/2 ; x0+winsize/2 ; x0+winsize/2 ; x0-winsize/2]; % Update x
        y = [y0-winsize/2 ; y0+winsize/2 ; y0+winsize/2 ; y0-winsize/2]; % Update y
        [XX,YY] = ndgrid([x(1):1:x(3)],[y(1):1:y(3)]);
        tempf = ImgRef([x(1):1:x(3)],[y(1):1:y(3)]);
        DfDx = Df.DfDx((x(1)-DfCropWidth):1:(x(3)-DfCropWidth), (y(1)-DfCropWidth):1:(y(3)-DfCropWidth));
        DfDy = Df.DfDy((x(1)-DfCropWidth):1:(x(3)-DfCropWidth), (y(1)-DfCropWidth):1:(y(3)-DfCropWidth));
    end

    meanf = mean(tempf(:)); bottomf = sqrt((length(tempf(:))-1)*var(tempf(:)));
    H2 = H(5:6,5:6)*2/(bottomf^2) + [mu 0; 0 mu];
    %%%%%%%%%%%%%%%% Tried below, not succeed %%%%%%%%%%%%%%%%%%%
    % if CrackOrNot == 0
    %     H2 = H(5:6,5:6)*2/(bottomf^2) + [mu 0; 0 mu];
    % else
    %     H2 = H*2/(bottomf^2) + [beta 0 0 0 0 0; 0 beta 0 0 0 0; 0 0 beta 0 0 0; 0 0 0 beta 0 0 ; 0 0 0 0 mu 0; 0 0 0 0 0 mu];
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % --------------------------
    % Initialize while loop
    normOfWOld=2; normOfWNew=1; normOfWNewAbs=1; stepwithinwhile=0;
    switch ICGNmethod   % For Gauss-Newton method
        case 'LevenbergMarquardt'
            delta = 0.001; % For Levenberg-Marquardt method
            KappaOld=1e10; KappaNew=1e10; KappaStore=zeros(10,1); PStore=zeros(10,6);
        otherwise % 'GaussNewton'
            delta = 0;
    end


    while( (stepwithinwhile<=100) && (normOfWNew>tol) && (normOfWNewAbs>tol) )

        stepwithinwhile = stepwithinwhile+1;

        % Find region for g
        tempCoordxMat = XX - x0*ones(winsize+1,winsize+1);
        tempCoordyMat = YY - y0*ones(winsize+1,winsize+1);
        u22 = (1+P(1))*tempCoordxMat + P(3)*tempCoordyMat + (x0+P(5))*ones(winsize+1,winsize+1);
        v22 = P(2)*tempCoordxMat + (1+P(4))*tempCoordyMat + (y0+P(6))*ones(winsize+1,winsize+1);

        row1 = find(u22<3); row2 = find(u22>imgSize(1)-2); row3 = find(v22<3); row4 = find(v22>imgSize(2)-2);
        if ~isempty([row1; row2; row3; row4])
            normOfWNew = 1e6;
            % warning('Out of image boundary!')
            break;
        else
            %tempg = imggNormalizedbc.eval(u22,v22)
            %tempg = ba_interp2(ImgDef, v22, u22, 'cubic');

            % BicubicBspline interpolation
            tempg = ba_interp2_spline(ImgDef, v22, u22, 'cubicspline');


            %%%%%%%%%% Visualize tempf and tempg %%%%%%%%%%
            %         stepwithinwhile
            %         subplot(1,2,1); surf(tempf,'edgecolor','none'); view(2); axis equal; axis tight;
            %         subplot(1,2,2); surf(tempg,'edgecolor','none'); view(2); axis equal; axis tight;
            %         pause;
            % ====== Old version codes ======
            % tempg = zeros(size(tempf,1)*size(tempf,2),1);
            % [tempCoordy, tempCoordx] = meshgrid(1:winsize+1,1:winsize+1);
            % tempCoordx = tempCoordx(:); tempCoordy = tempCoordy(:);
            %
            % for tempij = 1:size(tempCoordx,1)
            %     tempg(tempij,1)= ...
            %         fungInterpolation_g(  u22(tempCoordx(tempij),tempCoordy(tempij)) , v22(tempCoordx(tempij),tempCoordy(tempij)), ...
            %         g((floor(u22(tempCoordx(tempij),tempCoordy(tempij)))-1):(floor(u22(tempCoordx(tempij),tempCoordy(tempij)))+2), ...
            %         (floor(v22(tempCoordx(tempij),tempCoordy(tempij)))-1):(floor(v22(tempCoordx(tempij),tempCoordy(tempij)))+2))  );
            % end
            %
            % tempg = reshape(tempg, winsize+1, winsize+1);
            meang = mean(tempg(:));
            bottomg = sqrt((length(tempg(:))-1)*var(tempg(:)));

            % ============ For Levenberg-Marquardt method ============
            switch ICGNmethod
                case 'LevenbergMarquardt'
                    % Compute functinoal error
                    KappaOld = KappaNew;
                    Kappatemp = (tempf-meanf)/bottomf - (tempg-meang)/bottomg;
                    Kappatemp = Kappatemp.*Kappatemp;
                    KappaNew = sum(Kappatemp(:));

                    if KappaNew < 1.02*KappaOld
                        delta = delta/10;
                    else
                        delta = delta*10;
                        % Perform P inverse
                        DeltaP = -DeltaP;
                        tempP1 =  (-DeltaP(1)-DeltaP(1)*DeltaP(4)+DeltaP(2)*DeltaP(3))/tempDelP;
                        tempP2 =  -DeltaP(2)/tempDelP;
                        tempP3 =  -DeltaP(3)/tempDelP;
                        tempP4 =  (-DeltaP(4)-DeltaP(1)*DeltaP(4)+DeltaP(2)*DeltaP(3))/tempDelP;
                        tempP5 =  (-DeltaP(5)-DeltaP(4)*DeltaP(5)+DeltaP(3)*DeltaP(6))/tempDelP;
                        tempP6 =  (-DeltaP(6)-DeltaP(1)*DeltaP(6)+DeltaP(2)*DeltaP(5))/tempDelP;

                        tempMatrix = [1+P(1) P(3) P(5); P(2) 1+P(4) P(6); 0 0 1]*...
                            [1+tempP1 tempP3 tempP5; tempP2 1+tempP4 tempP6; 0 0 1];

                        P1 = tempMatrix(1,1)-1;
                        P2 = tempMatrix(2,1);
                        P3 = tempMatrix(1,2);
                        P4 = tempMatrix(2,2)-1;
                        P5 = tempMatrix(1,3);
                        P6 = tempMatrix(2,3);
                        P = [P1 P2 P3 P4 P5 P6]';
                    end

                    % Find region for g
                    % [tempCoordy, tempCoordx] = meshgrid(y(1):y(3),x(1):x(3));
                    %Repeated! tempCoordx = XX - x0*ones(winsize+1,winsize+1);
                    %Repeated! tempCoordy = YY - y0*ones(winsize+1,winsize+1);
                    u22 = (1+P(1))*tempCoordxMat + P(3)*tempCoordyMat + (x0+P(5))*ones(winsize+1,winsize+1);
                    v22 = P(2)*tempCoordxMat + (1+P(4))*tempCoordyMat + (y0+P(6))*ones(winsize+1,winsize+1);

                    tempg = ImgDef.eval(u22,v22);
                    % ====== Old version codes ======
                    % tempg = zeros(size(tempf,1)*size(tempf,2),1);
                    %
                    % [tempCoordy, tempCoordx] = meshgrid(1:winsize+1,1:winsize+1);
                    % tempCoordx = tempCoordx(:); tempCoordy = tempCoordy(:);
                    %
                    % parfor tempij = 1:size(tempCoordx,1)
                    %     tempg(tempij)= ...
                    %         fungInterpolation_g(u22(tempCoordx(tempij),tempCoordy(tempij)), v22(tempCoordx(tempij),tempCoordy(tempij)), ...
                    %         g(floor(u22(tempCoordx(tempij),tempCoordy(tempij)))-1:floor(u22(tempCoordx(tempij),tempCoordy(tempij)))+2, ...
                    %         floor(v22(tempCoordx(tempij),tempCoordy(tempij)))-1:floor(v22(tempCoordx(tempij),tempCoordy(tempij)))+2));
                    % end
                    %
                    % tempg = reshape(tempg, winsize+1, winsize+1);
                    % ==================================

                    meang = mean(tempg(:));
                    bottomg = sqrt((length(tempg(:))-1)*var(tempg(:)));

                otherwise
            end

            % ============ End of Levenberg-Marquardt method ============

            % Assemble b vector
            %         b = zeros(6,1);
            %         [tempCoordy, tempCoordx] = meshgrid(y(1):y(3),x(1):x(3));
            %         tempCoordx = tempCoordx(:); tempCoordy = tempCoordy(:);
            %
            %         % Compute for all the pixels in the master part
            %         for tempij = 1:size(tempCoordx,1)
            %             if PassCrackOrNot == 1 % Crack pass through this local subset
            %                 % if x0 < CrackTip(1) % Left part is master, and right part will be discarded
            %                 if CrackPathCen(1)*x0 + CrackPathCen(2)*y0 + 1 > 0  % Bottom part is master, and top part will be discarded
            %                     % if tempCoordx(tempij) > CrackTip(1) && tempCoordy(tempij) < CrackTip(2)
            %                     %     tempCoordx(tempij) = tempCoordx(tempij)-winsize;
            %                     % end
            %                     % if tempCoordx(tempij)+1-x(1)>1 && tempCoordx(tempij)>x0-0.8*winsize
            %                     %     b = b+ bottomf*([DfDx(tempCoordx(tempij)-DfCropWidth,tempCoordy(tempij)-DfCropWidth) DfDy(tempCoordx(tempij)-DfCropWidth,tempCoordy(tempij)-DfCropWidth)]*...
            %                     %         [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
            %                     %         ((tempf(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meanf)/bottomf - ...
            %                     %         (tempg(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meang)/bottomg);
            %                     % end
            %                     if tempCoordx(tempij)*CrackPathCen(1) + tempCoordy(tempij)*CrackPathCen(2) + 1 < 0
            %                         tempCoordy(tempij) = tempCoordy(tempij)-winsize;
            %                     end
            %                     if (tempCoordy(tempij)+1-y(1) > 1) && (tempCoordy(tempij)+1-y(1) < size(tempf,2)) ...
            %                         && (tempCoordx(tempij)*CrackPath2(1) + tempCoordy(tempij)*CrackPath2(2) + 1 > 0)
            %                            b = b+ bottomf*([DfDx(tempCoordx(tempij)-DfCropWidth,tempCoordy(tempij)-DfCropWidth) DfDy(tempCoordx(tempij)-DfCropWidth,tempCoordy(tempij)-DfCropWidth)]*...
            %                             [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
            %                             ((tempf(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meanf)/bottomf - ...
            %                             (tempg(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meang)/bottomg);
            %                     end
            %                 % elseif (x0 >= CrackTip(1)) % Right part is master, and left part will be discarded
            %                 elseif CrackPathCen(1)*x0 + CrackPathCen(2)*y0 + 1 <= 0  % Top part is master, and bottom part will be discarded
            %                     % if tempCoordx(tempij) < CrackTip(1) && tempCoordy(tempij) < CrackTip(2)
            %                     %     tempCoordx(tempij) = tempCoordx(tempij)+winsize;
            %                     % end
            %                     % if tempCoordx(tempij)<size(DfDx,1) && tempCoordx(tempij)+1-x(1)<size(tempf,1) && tempCoordx(tempij)<x0+0.5*winsize
            %                     %     b = b+ bottomf*([DfDx(tempCoordx(tempij)-DfCropWidth,tempCoordy(tempij)-DfCropWidth) DfDy(tempCoordx(tempij)-DfCropWidth,tempCoordy(tempij)-DfCropWidth)]*...
            %                     %         [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
            %                     %         ((tempf(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meanf)/bottomf - ...
            %                     %         (tempg(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meang)/bottomg);
            %                     % end
            %                     if tempCoordx(tempij)*CrackPathCen(1) + tempCoordy(tempij)*CrackPathCen(2) + 1 > 0
            %                         tempCoordy(tempij) = tempCoordy(tempij)+winsize;
            %                     end
            %                     if (tempCoordy(tempij)-DfCropWidth < size(DfDy,2)) && (tempCoordy(tempij)+1-y(1) < size(tempf,2)) ...
            %                         && (tempCoordx(tempij)*CrackPath1(1) + tempCoordy(tempij)*CrackPath1(2) + 1 < 0)
            %                             b = b + bottomf*([DfDx(tempCoordx(tempij)-DfCropWidth,tempCoordy(tempij)-DfCropWidth) DfDy(tempCoordx(tempij)-DfCropWidth,tempCoordy(tempij)-DfCropWidth)]*...
            %                             [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
            %                             ((tempf(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meanf)/bottomf - ...
            %                             (tempg(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meang)/bottomg);
            %                     end
            %                 end
            %             else
            %                 b = b + bottomf*([DfDx(tempCoordx(tempij)-DfCropWidth,tempCoordy(tempij)-DfCropWidth) DfDy(tempCoordx(tempij)-DfCropWidth,tempCoordy(tempij)-DfCropWidth)]*...
            %                     [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
            %                     ((tempf(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meanf)/bottomf - ...
            %                     (tempg(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meang)/bottomg);
            %             end
            %         end

            b2 = zeros(6,1);
            tempfMinustempg = (tempf-meanf*ones(winsize+1,winsize+1))/bottomf - (tempg-meang*ones(winsize+1,winsize+1))/bottomg;
            % b2(1) = sum(sum( (XX-x0).*DfDx.*tempfMinustempg ));
            % b2(2) = sum(sum( (XX-x0).*DfDy.*tempfMinustempg ));
            % b2(3) = sum(sum( (YY-y0).*DfDx.*tempfMinustempg ));
            % b2(4) = sum(sum( (YY-y0).*DfDy.*tempfMinustempg ));
            b2(5) = sum(sum( DfDx.*tempfMinustempg ));
            b2(6) = sum(sum( DfDy.*tempfMinustempg ));

            b = bottomf * b2;

            tempb = b(5:6)*2/(bottomf^2)  + [mu*(P(5)-UOld(1)-vdual(1)); mu*(P(6)-UOld(2)-vdual(2))];
            %%%%%%%%%%%%%%%% Tried below, not succeed %%%%%%%%%%%%%%%%%%%
            % if CrackOrNot == 0
            %     tempb = b(5:6)*2/(bottomf^2)  + [mu*(P(5)-UOld(1)-vdual(1)); mu*(P(6)-UOld(2)-vdual(2))];
            % else
            %     tempb = b*2/(bottomf^2) +[beta*(P(1)-FOld(1)-udual(1)); beta*(P(2)-FOld(2)-udual(2));
            %                               beta*(P(3)-FOld(3)-udual(3)); beta*(P(4)-FOld(4)-udual(4));
            %                               mu*(P(5)-UOld(1)-vdual(1)); mu *(P(6)-UOld(2)-vdual(2))];
            % end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            delta = 1e-3; DeltaP = [0 0 0 0 0 0];
            H2 = H(5:6,5:6)*2/(bottomf^2) + [mu 0; 0 mu];
            tempH = (H2 + delta*max(diag(H2))*eye(2));

            DeltaP(5:6) = -tempH\tempb;
            %%%%%%%%%%%%%%%% Tried below, not succeed %%%%%%%%%%%%%%%%%%%
            % if CrackOrNot == 0
            %     DeltaP(5:6) = -tempH\tempb;
            % else
            %     DeltaP = -(H2 + delta*diag(diag(H2))) \ tempb;
            % end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            normOfWOld = normOfWNew;
            normOfWNew = norm(tempb(:)); normOfWNewAbs = normOfWNew;

            if stepwithinwhile == 1
                normOfWNewInit = normOfWNew;
            end

            if normOfWNewInit > tol
                normOfWNew = normOfWNew/normOfWNewInit;
            else
                normOfWNew = 0;
            end

            if (normOfWNew<tol) || (normOfWNewAbs<tol)
                break
            else

                tempDelP =  ((1+DeltaP(1))*(1+DeltaP(4)) - DeltaP(2)*DeltaP(3));
                if (tempDelP ~= 0)
                    tempP1 =  (-DeltaP(1)-DeltaP(1)*DeltaP(4)+DeltaP(2)*DeltaP(3))/tempDelP;
                    tempP2 =  -DeltaP(2)/tempDelP;
                    tempP3 =  -DeltaP(3)/tempDelP;
                    tempP4 =  (-DeltaP(4)-DeltaP(1)*DeltaP(4)+DeltaP(2)*DeltaP(3))/tempDelP;
                    tempP5 =  (-DeltaP(5)-DeltaP(4)*DeltaP(5)+DeltaP(3)*DeltaP(6))/tempDelP;
                    tempP6 =  (-DeltaP(6)-DeltaP(1)*DeltaP(6)+DeltaP(2)*DeltaP(5))/tempDelP;

                    tempMatrix = [1+P(1) P(3) P(5); P(2) 1+P(4) P(6); 0 0 1]*...
                        [1+tempP1 tempP3 tempP5; tempP2 1+tempP4 tempP6; 0 0 1];

                    P1 = tempMatrix(1,1)-1;
                    P2 = tempMatrix(2,1);
                    P3 = tempMatrix(1,2);
                    P4 = tempMatrix(2,2)-1;
                    P5 = tempMatrix(1,3);
                    P6 = tempMatrix(2,3);
                    P = [P1 P2 P3 P4 P5 P6]';
                else
                    disp( 'Det(DeltaP)==0!' );
                    break
                end

            end
        end
    end % end of while
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (normOfWNew<tol) || (normOfWNewAbs<tol)
    else
        stepwithinwhile = 101;
    end

    if (isnan(normOfWNew)==1)
        stepwithinwhile = 101;
    end


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else

    stepwithinwhile = 102;
    H = zeros(6,6);

end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U = [P(5);P(6)]; H = H;
% F(1) = P(1); F(2) = P(2); F(3) = P(3); F(4) = P(4);

if nargout >= 4
    varargout{1} = BicubicBspline;
end

end


