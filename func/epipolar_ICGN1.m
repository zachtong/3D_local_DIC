function  ZNSSD = epipolar_ICGN1(pointCandidate,Df,ImgRef,ImgDef,winsize,tol,ICGNmethod)
% Special ICGN-1st_order version for epipolar searching. (Fix U,V)
% ----------------------------------------------
%   INPUT:
%          pointCandidate
%          Df                   Image grayscale value gradients
%          ImgRef               Reference image
%          ImgDef               Deformed image
%          winsize              DIC parameter subset size
%          ICGNmethod           ICGN iteration scheme: 'GaussNewton' -or- 'LevenbergMarquardt'
%          tol                  ICGN iteration stopping threshold
%
%   OUTPUT: U                   Disp vector: [Ux_node1, Uy_node1, ... , Ux_nodeN, Uy_nodeN]';
%           F                   Deformation gradient tensor
%                               F = [F11_node1, F21_node1, F12_node1, F22_node1, ... , F11_nodeN, F21_nodeN, F12_nodeN, F22_nodeN]';
%           ZNSSD
% ==============================================

%% Initialization
warning('off');
DfDxStartx = Df.DfAxis(1); DfDxStarty = Df.DfAxis(3);
imgSize = Df.imgSize;
x0 = pointCandidate(1);
y0 = pointCandidate(2);

%% ---------------------------
% Find local subset region
x = [x0-winsize/2 ; x0+winsize/2 ; x0+winsize/2 ; x0-winsize/2];  % [coordinates(elements(j,:),1)];
y = [y0-winsize/2 ; y0+winsize/2 ; y0+winsize/2 ; y0-winsize/2];  % [coordinates(elements(j,:),2)];

% ---------------------------
% Initialization: Get P0
P0 = [0 0 0 0 0 0]'; % P0(5) P0(6) is fixed
P = P0;

% ---------------------------
% Find region for f
[XX,YY] = ndgrid(x(1):1:x(3),y(1):1:y(3));
%tempf = imgfNormalizedbc.eval(XX,YY);
%DfDx = imgfNormalizedbc.eval_Dx(XX,YY);
%DfDy = imgfNormalizedbc.eval_Dy(XX,YY);
tempf = ImgRef(x(1):1:x(3),y(1):1:y(3));
DfDx = Df.DfDx((x(1)-DfDxStartx):1:(x(3)-DfDxStartx), (y(1)-DfDxStarty):1:(y(3)-DfDxStarty));
DfDy = Df.DfDy((x(1)-DfDxStartx):1:(x(3)-DfDxStartx), (y(1)-DfDxStarty):1:(y(3)-DfDxStarty));


%% %%%%%%%% If there are >50% of the subset are painted with patterns %%%%%%%%%%%%
[DfDxImgMaskIndRow,~] = find(DfDx==0);
if length(DfDxImgMaskIndRow)<0.50*(winsize+1)^2

    if length(DfDxImgMaskIndRow)>0.1*(winsize+1)^2
        winsize = max(round(sqrt(2)*winsize));
        x = [x0-winsize/2 ; x0+winsize/2 ; x0+winsize/2 ; x0-winsize/2]; % Update x
        y = [y0-winsize/2 ; y0+winsize/2 ; y0+winsize/2 ; y0-winsize/2]; % Update y
        [XX,YY] = ndgrid([x(1):1:x(3)],[y(1):1:y(3)]);
        tempf = ImgRef([x(1):1:x(3)],[y(1):1:y(3)]);
        DfDx = Df.DfDx((x(1)-DfDxStartx):1:(x(3)-DfDxStartx), (y(1)-DfDxStarty):1:(y(3)-DfDxStarty));
        DfDy = Df.DfDy((x(1)-DfDxStartx):1:(x(3)-DfDxStartx), (y(1)-DfDxStarty):1:(y(3)-DfDxStarty));
    end

    H2 = zeros(6,6); DfDxSq = (DfDx.^2); DfDySq = (DfDy.^2); DfDxDfDy = DfDx.*DfDy;
    XXSq = (XX-x0).^2; YYSq = (YY-y0).^2; XXYY = (XX-x0).*(YY-y0);
    H2(1,1) = sum(sum(XXSq.*DfDxSq));       H2(1,2) = sum(sum(XXSq.*DfDxDfDy ));
    H2(1,3) = sum(sum( XXYY.*DfDxSq ));     H2(1,4) = sum(sum( XXYY.*DfDxDfDy ));
    H2(1,5) = sum(sum( (XX-x0).*DfDxSq ));  H2(1,6) = sum(sum( (XX-x0).*DfDxDfDy ));
    H2(2,2) = sum(sum(XXSq.*DfDySq));       H2(2,3) = H2(1,4);
    H2(2,4) = sum(sum( XXYY.*DfDySq ));     H2(2,5) = H2(1,6);
    H2(2,6) = sum(sum( (XX-x0).*DfDySq ));  H2(3,3) = sum(sum( YYSq.*DfDxSq ));
    H2(3,4) = sum(sum( YYSq.*DfDxDfDy ));   H2(3,5) = sum(sum( (YY-y0).*DfDxSq ));
    H2(3,6) = sum(sum( (YY-y0).*DfDxDfDy ));H2(4,4) = sum(sum( YYSq.*DfDySq ));
    H2(4,5) = H2(3,6);  H2(4,6) = sum(sum((YY-y0).*DfDySq)); H2(5,5) = sum(sum(DfDxSq));
    H2(5,6) = sum(sum(DfDxDfDy)); H2(6,6) = sum(sum(DfDySq));
    H = H2 + H2' - diag(diag(H2));

    meanf = mean(tempf(:));
    bottomf = sqrt((length(tempf(:))-1)*var(tempf(:))); % Sum of std_deviation

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

    while( (stepwithinwhile<=30) && (normOfWNew>tol) && (normOfWNewAbs>tol) ) % Iterate 30 times

        stepwithinwhile = stepwithinwhile+1;

        % Find region for g
        % %[tempCoordy, tempCoordx] = meshgrid(y(1):y(3),x(1):x(3));
        tempCoordxMat = XX - x0*ones(winsize+1,winsize+1);
        tempCoordyMat = YY - y0*ones(winsize+1,winsize+1);
        u22 = (1+P(1))*tempCoordxMat + P(3)*tempCoordyMat + (x0+P(5))*ones(winsize+1,winsize+1);
        v22 = P(2)*tempCoordxMat + (1+P(4))*tempCoordyMat + (y0+P(6))*ones(winsize+1,winsize+1);

        row1 = find(u22<3); row2 = find(u22>imgSize(1)-2); row3 = find(v22<3); row4 = find(v22>imgSize(2)-2);
        if ~isempty([row1; row2; row3; row4])
            normOfWNew = 1e6; % warning('Out of image boundary!')
            break;
        else

            %tempg = imggNormalizedbc.eval(u22,v22)
            tempg = ba_interp2(ImgDef, v22, u22, 'cubic');

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
                        tempP1 =  (-DeltaP(1)-DeltaP(1)*DeltaP(4)+DeltaP(2)*DeltaP(3))/temp;
                        tempP2 =  -DeltaP(2)/temp;
                        tempP3 =  -DeltaP(3)/temp;
                        tempP4 =  (-DeltaP(4)-DeltaP(1)*DeltaP(4)+DeltaP(2)*DeltaP(3))/temp;
                        tempP5 =  (-DeltaP(5)-DeltaP(4)*DeltaP(5)+DeltaP(3)*DeltaP(6))/temp;
                        tempP6 =  (-DeltaP(6)-DeltaP(1)*DeltaP(6)+DeltaP(2)*DeltaP(5))/temp;

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
                    meang = mean(tempg(:));
                    bottomg = sqrt((length(tempg(:))-1)*var(tempg(:)));

                otherwise
            end

            % ====== Assemble b vector fast version ======
            b2 = zeros(6,1);
            tempfMinustempg = (tempf-meanf*ones(winsize+1,winsize+1))/bottomf - (tempg-meang*ones(winsize+1,winsize+1))/bottomg;
            b2(1) = sum(sum( (XX-x0).*DfDx.*tempfMinustempg ));
            b2(2) = sum(sum( (XX-x0).*DfDy.*tempfMinustempg ));
            b2(3) = sum(sum( (YY-y0).*DfDx.*tempfMinustempg ));
            b2(4) = sum(sum( (YY-y0).*DfDy.*tempfMinustempg ));
            b2(5) = sum(sum( DfDx.*tempfMinustempg ));
            b2(6) = sum(sum( DfDy.*tempfMinustempg ));

            b = bottomf * b2;

            normOfWOld = normOfWNew;
            normOfWNew = norm(b(:)); normOfWNewAbs = normOfWNew;

            if stepwithinwhile ==1
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
                % DeltaP = [0 0 0 0 0 0];
                % tempH = (H + delta*diag(diag(H)));
                % tempb = b;
                % DeltaP(5:6) = -tempH(5:6,5:6)\tempb(5:6);
                DeltaP = -(H + delta*diag(diag(H))) \ b;
                detDeltaP =  ((1+DeltaP(1))*(1+DeltaP(4)) - DeltaP(2)*DeltaP(3));
                if (detDeltaP ~= 0)
                    tempP1 =  (-DeltaP(1)-DeltaP(1)*DeltaP(4)+DeltaP(2)*DeltaP(3))/detDeltaP;
                    tempP2 =  -DeltaP(2)/detDeltaP;
                    tempP3 =  -DeltaP(3)/detDeltaP;
                    tempP4 =  (-DeltaP(4)-DeltaP(1)*DeltaP(4)+DeltaP(2)*DeltaP(3))/detDeltaP;
                    tempP5 =  (-DeltaP(5)-DeltaP(4)*DeltaP(5)+DeltaP(3)*DeltaP(6))/detDeltaP;
                    tempP6 =  (-DeltaP(6)-DeltaP(1)*DeltaP(6)+DeltaP(2)*DeltaP(5))/detDeltaP;

                    tempMatrix = [1+P(1) P(3) P(5); P(2) 1+P(4) P(6); 0 0 1]*...
                        [1+tempP1 tempP3 tempP5; tempP2 1+tempP4 tempP6; 0 0 1];

                    P1 = tempMatrix(1,1)-1;
                    P2 = tempMatrix(2,1);
                    P3 = tempMatrix(1,2);
                    P4 = tempMatrix(2,2)-1;
                    P5 = 0;
                    P6 = 0; % P5 P6 fixed
                    P = [P1 P2 P3 P4 P5 P6]';
                else
                    disp(['Det(DeltaP)==0!'])
                    break
                end

            end
        end
    end % end of while
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (normOfWNew<tol) || (normOfWNewAbs<tol)
        % elementsLocalMethodConvergeOrNot = 1;
    else
        stepwithinwhile = 101;
    end

    if (isnan(normOfWNew)==1)
        stepwithinwhile = 101;
    end


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if norm(diag(H)) > abs(eps)

    H = zeros(6,6);
    stepwithinwhile = 102;

end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U(1) = P(5); U(2) = P(6);
F(1) = P(1); F(2) = P(2); F(3) = P(3); F(4) = P(4);
% HGlobal = [H(1:6) H(8:12) H(15:18) H(22:24) H(29:30) H(36)];

% ZNSSD
Kappatemp = (tempf-meanf)/bottomf - (tempg-meang)/bottomg;
Kappatemp = Kappatemp.*Kappatemp;
KappaNew = sum(Kappatemp(:));
ZNSSD = KappaNew;



end

