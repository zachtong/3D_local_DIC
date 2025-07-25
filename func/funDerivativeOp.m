function A = funDerivativeOp(M,N,h)
% FUNCTION A = funDerivativeOp(M,N,h)
% To generate a first order gradient operator {A} in the 2D case such that 
% {F} = {A}{U} where displacement vector U = [Ux_node1, Uy_node1, ... , Ux_nodeN, Uy_nodeN]',
% and deformation gradient F = [F11_node1, F21_node1, F12_node1, F22_node1, ... 
% ... , F11_nodeN, F21_nodeN, F12_nodeN, F22_nodeN]';
% ----------------------------------------------
%
%   INPUT: M,N   Mesh size in the x- and y-directions
%
%   OUTPUT: A    Finite difference operator
%
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 2018.03, 2020.12 
% ============================================== 


%% Initialization
DIM = 2;

A = sparse(DIM^2*M*N, DIM*M*N); 
sizeA = size(A); 

[XX,YY] = ndgrid(1:M,1:N); XX = XX(:); YY = YY(:);  
INDEXI = ones(2*DIM^2*M*N ,1); INDEXJ = INDEXI; INDEXVAL = INDEXI;


%% ------ Inside assign values ------
for tempi = 1:M*N %%%%% Write into one for-loop in MATLAB instead of two for-loops
    tempx = XX(tempi); tempy = YY(tempi);  

            % Determine BorderOrNot
            BorderOrNot = ones(1,DIM*2);  
            if tempy == 1, BorderOrNot(2) = 0; end; if tempy == N, BorderOrNot(1) = 0; end
            if tempx == 1, BorderOrNot(4) = 0; end; if tempx == M, BorderOrNot(3) = 0; end
            
            index = tempx + M*(tempy-1)  ; % Find the index of the point (tempx,tempy,tempz);
            indexNeighbors = index*ones(1,DIM*2)+BorderOrNot.*[M,-M,1,-1]; %{Back,Front,Right,Left,Top,Bottom}
            % indexBack = index+M; indexFront = index-M; indexLeft = index-1;
            % indexRight = index+1; indexTop = index+M*N; indexBottom = index-M*N;
            
            % Find index of affine deformation gradient tensor {F}
            indexFrow = DIM^2*index*ones(1,DIM^2)+[-(DIM^2-1):1:0]; %{F11,F21,F31,F12,F22,F32,F13,F23,F33};
            
            
            indexFcol1 = [DIM*indexNeighbors(4)*ones(1,DIM)+[-(DIM-1):1:0], DIM*indexNeighbors(2)*ones(1,DIM)+[-(DIM-1):1:0] ];
            indexFcol2 = [DIM*indexNeighbors(3)*ones(1,DIM)+[-(DIM-1):1:0], DIM*indexNeighbors(1)*ones(1,DIM)+[-(DIM-1):1:0] ];
   

            
            INDEXI(2*4*tempi-7:2*4*tempi-0) = [indexFrow,indexFrow]';
            INDEXJ(2*4*tempi-7:2*4*tempi-0) = [indexFcol1,indexFcol2]';
            INDEXVAL(2*4*tempi-7:2*4*tempi-4) = -ones(4,1);
            
            if BorderOrNot(3)*BorderOrNot(4) == 0, INDEXVAL(2*4*tempi-7:2*4*tempi-6) = -2*ones(2,1); INDEXVAL(2*4*tempi-3:2*4*tempi-2) = 2*ones(2,1); end
            if BorderOrNot(1)*BorderOrNot(2) == 0, INDEXVAL(2*4*tempi-5:2*4*tempi-4) = -2*ones(2,1); INDEXVAL(2*4*tempi-1:2*4*tempi ) = 2*ones(2,1); end
             
            
 end

A = (1.0/(2*h)) * sparse(INDEXI,INDEXJ,INDEXVAL, sizeA(1),sizeA(2));
end