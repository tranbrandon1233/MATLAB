%Script: Final Project
%Description:  The goal of this project is to find the unknown temperature distribution of an
% aluminum fin that is subjected to a heat input using the point
% collocation method and the radial basis function as the interpolant
%Author: Brandon Tran
%UID: 705830462

clc;clearvars; close all;  %Clear vars from previous programs
rng('shuffle');

%Initialize vars
eps = 1;
tol = 0.3;
nPlot = 9999;
%Vars for 1D
L1D = 10;
nNodes = 10;
%Vars for 2D
L2D = 16;
H = 9;
k=190;
Tbar = @(x,y) x*y/32;
nNodesXOmega = 14;
nNodesYOmega = 7;
nNodesXGamma = 8;
nNodesYGamma = 15;

% Define nodes on Gamma for 2D problems:
nNodesGamma = 2*(nNodesYGamma+nNodesXGamma);
%Initialize x- and y-values of nodes on Gamma
nodesGamma = zeros(nNodesGamma,2);

%Node distance in the x- and y-direction (on Gamma)
dxGamma = L2D/nNodesXGamma;
dyGamma = H/nNodesYGamma;

%Place nodes on boundary Gamma
% Bottom part: x-values of nodes on Gamma
    xGammaBottom = linspace(0,L2D-dxGamma,nNodesXGamma);
% Right part: y-values of nodes on Gamma
    yGammaRight = linspace(0,H-dyGamma,nNodesYGamma);
% Top part: x-values of nodes on Gamma
    xGammaTop = linspace(0,L2D-dxGamma,nNodesXGamma);
% Left part: y-values of nodes on Gamma
    yGammaLeft = linspace(0,H-dyGamma,nNodesYGamma);

    %Assign x and y coordinates to nodes
    %Bottom: Nodes on Gamma [x-coord. = xGammaBottom, y-coord. = 0]
    nodesGamma(1:nNodesXGamma,1) = xGammaBottom;
    nodesGamma(1:nNodesXGamma,2) = zeros(nNodesXGamma,1);
    %Right part: Nodes on Gamma [x-coord. = L2D, y-coord. = yGammaRight]
    nodesGamma(nNodesXGamma+1:nNodesXGamma+nNodesYGamma,1) = L2D*ones(nNodesYGamma, 1);
    nodesGamma(nNodesXGamma+1: nNodesXGamma + nNodesYGamma, 2) = yGammaRight;
    %Top part:  Nodes on Gamma [x-coord. = H, y-coord. = xGammaTop]
    nodesGamma(nNodesXGamma+nNodesYGamma+1:2*nNodesXGamma+nNodesYGamma,1) = H*ones(nNodesXGamma, 1);
    nodesGamma(nNodesXGamma+nNodesYGamma+1:2*nNodesXGamma + nNodesYGamma, 2) = xGammaTop;
    %Left: Nodes on Gamma [x-coord. = yGammaLeft, y-coord. = 0]
    nodesGamma(2*nNodesXGamma+nNodesYGamma+1:end,1) = yGammaLeft;
    nodesGamma(2*nNodesXGamma+nNodesYGamma+1:end,2) = zeros(nNodesYGamma,1);
 
   %Ask user how nodes should be distributed
   randomInput = input("Should the nodes be equally (1) or randomly (2) distributed? (Enter 1 or 2): ");
   switch(randomInput)  %Set random flag
       case 1
           random = false;
       case 2
           random = true;
       otherwise
           error("Invalid input.")
   end

   %Distribute nodes
   [~,nodes]=defineNodes(nNodesXGamma,nNodesYGamma,tol,L1D,H,1,random);
    
   %Define anonymous function
   f=@(x) exp(x/8.*sin(2.*x-1));
   %Solve for array w
    w=solveLinearSystem(nNodes,nNodes,nodes,nodes,f,eps,k,1);
    %Plot results
    plotResults(nNodes,nPlot,nodes,w,f,eps,L1D,H,1)

    %Repeat for part 2

   %Ask user how nodes should be distributed
   randomInput = input("Should the nodes be equally (1) or randomly (2) distributed? (Enter 1 or 2): ");
   switch(randomInput)  %Set random flag
       case 1
           random = false;
       case 2
           random = true;
       otherwise
           error("Invalid input.")
   end

      %Distribute nodes
   [nNodesOmega,nodesOmega]=defineNodes(nNodesXOmega,nNodesYOmega,tol,L2D,H,2,random);
    
   %Define anonymous function
   f=Tbar;
   %Solve for array w
    w=solveLinearSystem(nNodesOmega,nNodesGamma,nodesOmega,nodesGamma,f,eps,k,2);
    %Plot results with appropriate parameters
    nNodes = nNodesYOmega*nNodesXOmega+nNodesGamma;
    nodes = [nodesOmega;nodesGamma];
    plotResults(nNodes,nPlot,nodes,w,f,eps,L2D,H,2)

   %Repeat for part 3

   %Ask user how nodes should be distributed
   randomInput = input("Should the nodes be equally (1) or randomly (2) distributed? (Enter 1 or 2): ");
   switch(randomInput)  %Set random flag
       case 1
           random = false;
       case 2
           random = true;
       otherwise
           error("Invalid input.")
   end

      %Distribute nodes
   [nNodesOmega,nodesOmega]=defineNodes(nNodesXOmega,nNodesYOmega,tol,L2D,H,3,random);
    
   %Define anonymous function
   f=Tbar;
   %Solve for array w
    w=solveLinearSystem(nNodesOmega,nNodesGamma,nodesOmega,nodesGamma,f,eps,k,3);
    %Plot results with appropriate parameters
    nNodes = nNodesYOmega*nNodesXOmega+nNodesGamma;
    nodes = [nodesOmega;nodesGamma];
    plotResults(nNodes,nPlot,nodes,w,f,eps,L2D,H,2)

function [nNodes,nodes] = defineNodes(nNodesX,nNodesY,tol,L,H, whichPart,random)
%Returns number and locatins of nodes based on the dimension determined by
%the number of nodes in the x and y directions, the tolerance tol of a
%random distribution of ndoes, the length L and height H, the part of the
%problem determined by whichPart, and the random flag to determine if
%random node generation should be used
    nNodes=0;
    if whichPart == 1
        %Check if function should distribute nodes equally or randomly
        if ~random
            nodes = linspace(0,L);
        else
            nodes =tol + (L - 2*tol) * rand(1, nNodesX); %randomly distributed nodes
            nodes(1) = 0;  %Ensure first value is 0 and last value is always L
            nodes(length(nodes)) = L;
        end
        %Handle random distribution
        nodes(2:length(nodes)-1) = mergeSort(nodes(2:length(nodes)-1));

        %Initialize parameters to loop through nodes
        nodeMoved = true;
        iterCount = 0;

        while nodeMoved
            nodeMoved = false;
            for i=2:nNodesX-2
                if abs(nodes(i+1)-nodes(i))<tol
                    nodes(i) = tol + (L - 2*tol) * rand(1); %Randomly redefine nodes
                    nodes(2:length(nodes)-1) = mergeSort(nodes(2:length(nodes)-1)); %Sort nodes again
                    nodeMoved = true;
                    iterCount = iterCount+1;
                end
            end
        end  
    else  %When whichPart is not 1
        nNodes = nodesX*nNodesY;  %Initialize vars
        nodes = zeros(nNodes*2,nNodes*2);
        if ~random  %If nodes are not random
            dxOmega = L/(nNodesX+1); %Use length for distribution based on x
            dyOmega = H/(nNodesY+1); %Use height for distribution based on y
            xOmega = linspace(dxOmega,L-dxOmega,nNodesX); %x of nodes in omega
            yOmega = linspace(dyOmega,H-dyOmega,nNodesY); %y of nodes in omega
            for i=1:length(xOmega)
                for j=1:length(yOmega)
                    row = (i-1)*nNodesY+j;
                    nodes(row,1) = xOmega(i);
                    nodes(row,2) = yOmega(j);
                end
            end
        else   %If randomly distributed nodes
            for i=1:nNodes %Loop over all nodes
                nodes(i,1) = tol+(L-2*tol)*rand;   % random x-values
                nodes(i,2) = tol+(H-2*tol)*rand;   % random y-values
            end
        end
    end
    %Initialize vars
    nodeMoved = true;
    iterCount = 0;
    while nodeMoved
        nodeMoved = false;
        for i = 1:nNodes
            for j=1:nNodes
                if i ~= j
                    dist = sqrt((nodes(i,1)-nodes(j,1))^2+(nodes(i,2)-nodes(j,2))^2);
                    % Check if node distance is less than tolerance
                    if dist < tol
                        %Create new values for nodes
                        nodes(i,1) = tol+(L-2*tol)*rand;
                        nodes(i,2) = tol+(H-2*tol)*rand;
                        nodeMoved = true;
                        iterCount = iterCount+1;
                        break;
                    end
                end
            end
        end
    end
    disp(iterCount)
end

function sol = solveLinearSystem(nNodesOmega,nNodesGamma,nodesOmega,nodesGamma,f,eps,k,whichPart)
%Finds the vector of unknown weights from the system A sol = b based on the
%system matrix A created using the basis functions (which use eps), the
%nodes from nodesOmega and nodesGamma, the heat conductivity k if whichPart
%is 3, and finding the vector b using f and Tbar if whichPart is 3
    if whichPart == 1 %Check whichPart
        %Initialize vars
        A = size(nNodesOmega,nNodesOmega);
        for i=1:length(A)
            for j=1:min(size(A))
                A(i,j) = evaluatePhiJatI(nodesOmega(j,:), nodesOmega(i,:), eps, 0);  %Find matrix A
            end
        end
        b = f(nodesOmega);  %Find b vector
    else %In all other cases,
        %Create submatrices of A
         AOmegaOmega = zeros(nNodesOmega,nNodesOmega);  
         AOmegaGamma = zeros(nNodesOmega,nNodesGamma);  
         AGammaOmega = zeros(nNodesGamma,nNodesOmega);
         AGammaGamma = zeros(nNodesGamma,nNodesGamma);
         %Find deriv to use
         if whichPart == 2
             deriv = 0;
         elseif whichPart == 3
             deriv = 2;
         end
         %Fill the submatrices
        %AOmegaOmega
        for i=1:nNodesOmega
            for j=1:nNodesOmega
                AOmegaOmega(i,j) = evaluatePhiJatI(nodesOmega(j,:), nodesOmega(i,:), eps, deriv);
            end
        end
        %AOmegaGamma
        for i=1:nNodesOmega
            for j=1:nNodesGamma
                AOmegaGamma(i,j) = evaluatePhiJatI(nodesGamma(j,:), nodesOmega(i,:), eps, deriv);
            end
        end
        %AGammaOmega
        for i=1:nNodesGamma
            for j=1:nNodesOmega
                AGammaOmega(i,j) = evaluatePhiJatI(nodesOmega(j,:), nodesGamma(i,:), eps, deriv);
            end
        end
        %AGammaGamma
        for i=1:nNodesGamma
            for j=1:nNodesGamma
                AGammaGamma(i,j) = evaluatePhiJatI(nodesGamma(j,:), nodesGamma(i,:), eps, deriv);
            end
        end
        %Find A matrix depending on whichPart parameter
        if whichPart == 2  %Part 2 definition for A
            A = [AOmegaOmega AOmegaGamma; AGammaOmega AGammaGamma];
            b = [f(nNodesOmega(:, 1));f(nNodesGamma(:, 1))];
        else  %Part 3 definition for A
             A = [-k*AOmegaOmega -k*AOmegaGamma; AGammaOmega AGammaGamma];
            b = [1000*f(nNodesOmega(:, 1));Tbar(nNodesGamma(:, 1))];
        end
    end
    sol = A\b; %Solve for linear system of equations
end

function phiVal = evaluatePhiJatI(nodeJ,nodeI,eps,deriv)
%Evaluates phi_j, the gradient if deriv = 1, or the divergence of the
%gradient for deriv = 2 based on the shape parameter eps evaluated at node
%x_i (or using any other node position)
    %Initialize vars based on input
    dim = size(nodeI,2);
    r=0;

    for i=1:dim
        r=r+(nodeI(i)-nodeJ(i))^2;  %Sum all squared differences to r
    end
    r  = sqrt(r); 
    if deriv == 0
        phiVal = sqrt(1+(eps*r)^2);
    elseif deriv == 1
        phiVal = 2*eps*r/(2*sqrt(1+(eps*r)^2));
    elseif deriv == 2
        phiVal = -(exp(-((y*x)/64-sin(x/2-8))^2)*(((98304*y^2*cos(x/2-8)-3072*y^3)*x+(196608*y^2-6291456*y*cos(x/2-8))*sin(x/2-8))*sin((y*x)/12)+((-9216*y^2*cos(x/2-8)^2+576*y^3*cos(x/2-8)-9*y^4)*x^2+(1179648*y*cos(x/2-8)^2-73728*y^2*cos(x/2-8)+1152*y^3+294912*y)*sin(x/2-8)*x+(-37748736*cos(x/2-8)^2+2359296*y*cos(x/2-8)-36864*y^2-18874368)*sin(x/2-8)^2+18874368*cos(x/2-8)^2-1179648*y*cos(x/2-8)+280576*y^2)*cos((y*x)/12)))/37748736 + (x^2*exp(-((x*y)/64-sin(x/2-8))^2)*((3072*x*y-196608*sin(x/2-8))*sin((x*y)/12)+(9*x^2*y^2-1152*sin(x/2-8)*x*y+36864*sin(x/2-8)^2-280576)*cos((x*y)/12)))/37748736;
    end
end

function plotResults(nNodes,nPlot,nodes,w,f,eps,L,H,whichPart)
%Creates figures with plots of the results using nNodes as the total number
%of nodes, nPlot points in each direction, the nodes array, the weights w,
%the function f, the shape parameter eps, the length L, the height H, and
%whichPart as the part number

    %Initialize variables
    xPlot = linspace(0,L,nPlot); % Define nodes (1D), x-value of nodes (2D)
    f_hAtNodes = zeros(nNodes,1);

    if whichPart == 1  %Part 1
        f_hAtXPlot = zeros(nPlot,1);   %Vars to track f_h and phi_j values in xPlot
        phi_jAtXPlot = zeros(nPlot,nNodes);
        
        for j=1:nNodes
            for i=1:nPlot
                f_hAtXPlot(i) = f_hAtXPlot(i)+w(j)*evaluatePhiJatI(nodes(j,:), nodes(i,:), eps, 1);
                phi_jAtXPlot(i,j) = phi_jAtXPlot(i,j)+evaluatePhiJatI(nodes(j,:), nodes(i,:), eps, 1);
            end

            for i=1:nNodes
                f_hAtNodes(i) = f_hAtXPlot(i); %Determine f_h-values at the nodes xi (i.e. f_h(xi))
            end
        end
        %Plot figures
        
        figure
        plot(f_hAtXPlot, xPlot)
        plot(f(xPlot), xPlot)
        scatter(nodes(:,1),f_hAtNodes(:,1))

        figure
        errorXPlot = f(xPlot)-f_hAtXPlot;
        errorXi = f(nodes) - f_hAtNodes(nodes);
        plot(errorXPlot,xPlot)
        plot(errorXi, nodes)
        plot(phi_jAtXPlot, xPlot)
    else  %Plots for parts 2 and 3
        yPlot = linspace(0,H,nPlot); %y value for nodes
        [X,Y] = meshgrid(xPlot,yPlot);
        fAtYX = f(X,Y);
        f_hAtYPlotXPlot = zeros(nPlot,nPlot);

        for j=1:nNodes
            for k=1:nPlot
                for l=1:nPlot
                    %Find f_h vales at points (yPlot,xPlot)
                    f_hAtYPlotXPlot(k,l) = f_hAtYPlotXPlot(k,l) + w(j)*evaluatePhiJatI(nodes(j,:),[xPlot(l),yPlot(k)], eps, 1);
                end
            end
        end
        for i=1:nNodes
                %Find f_h values at nodes (x_i,y_i)
                f_hAtNodes(i) = f_hAtYPlotXPlot(nodes(:,1),nodes(:,2));
        end
    end
    if whichPart == 3 %Check part for heat conduction problem
        %Create heat flux vector
        q_h=zeros(nNodes,2);
        for j=1:nNodes
            for i=1:nodes
                q_h(i,:) = q_h(i,:)+w(j)*evaluatePhiJatI(nodes(j,:), nodes(i,:), eps, 1);      
            end
        end
        q_h=-k*q_h; %Fourier's law
    end
    %Plot figures
    figure
    %Contour plot of f_h(yPlot,xPlot)
    contour(X,Y,f_hAtYPlotXPlot)
    %Plot nodes (x_i,y_i)
    scatter(nodes(:,1),nodes(:,2)) 
    if whichPart == 3
%          Plot vectorial q_h-values at nodes (x_i,y_i) if part 3
          quiver(nodes(:,1),nodes(:,2),q_h(:,1),q_h(:,2));
    end
    figure                            
%        Surface plot of f_h(yPlot,xPlot) (T_h(yPlot,xPlot) in Part 3)
          surf(X,Y,f_hAtYPlotXPlot)
%        Plot data points (x_i,y_i,f_h(x_i,y_i)) (T_h(x_i,y_i) in Part 3)
          scatter3(nodes(:,1),nodes(:,2),f_hAtNodes);
%       Other plots for Part 2
         if whichPart == 2
            figure                         
%         Contour plot of error function
            error = fAtYX-f_hAtYPlotXPlot;
            contour(X,Y,error)
%         Plot nodes (x_i,y_i)
           scatter(nodes(:,1),nodes(:,2))

            figure                          
%         Surface plot of error function e(yPlot,xPlot)
            surf(X,Y,error)
%         - Plot data points (x_i,y_i,e(x_i,y_i))
            scatter3(nodes(:,1),nodes(:,2),error(nodes(:,1),nodes(:,2)))
         end
end

% Merge sort functions from lecture
%====================================================================   
%%  Function - Merge Sort
function y = mergeSort(x)
%   Given a 1-D array x, return the array y with
%   entries in x sorted from smallest to largest
%   Initialization
n = length(x);
%   Error checks
if n < 1   % Input vector with less than 1 entry
    error('Input vector requires more entries!');
end
for i = 1:n
    if isnumeric(x(i)) == 0 % Non-numeric input vector
        error('Non-numeric entry in input vector!');
    end
end
if n == 1   % Termination condition
    y = x;
else
    %   Split input vector
    m   = floor(n/2);
    xL  = x(1:m);
    xR  = x(m+1:n);
    %   Sort each half
    yL  = mergeSort(xL);
    yR  = mergeSort(xR);
    %   Merge the sorted vectors
    y   = merge(yL,yR);
end
end

%====================================================================  
%%  Function - Merging Two Sorted Vectors
function z = merge(x,y)
%   Given two sorted 1-D arrays x and y, return the array z with
%   entries sorted from smallest to largest
%   Initialization
nx = length(x);
ny = length(y);
z = zeros(nx+ny,1);     % Initialize z with size nx + ny
% Initialize x,y pointer at left and initialize z counter
ixptr = 1;
iyptr = 1;
zcnt = 1;
%   Error checks
if nx < 1 || ny < 1  % Input vector with less than 1 entry
    error('Input vector requires more entries!');
end
for i = 1:nx
    if isnumeric(x(i)) == 0 % Non-numeric input vector x
        error('Non-numeric entry in input vector x!');
    end
end
for j = 1:ny
    if isnumeric(y(j)) == 0 % Non-numeric input vector y
        error('Non-numeric entry in input vector y!');
    end
end
%   Merge x and y into z vector by choosing the smallest value at each step
while ixptr <= nx && iyptr <= ny
    if x(ixptr) <= y(iyptr)
        z(zcnt) = x(ixptr);
        zcnt    = zcnt + 1;
        ixptr   = ixptr + 1;
    else
        z(zcnt) = y(iyptr);
        zcnt    = zcnt + 1;
        iyptr   = iyptr + 1;
    end
end
%   Deal with remaining entries in x or y
while ixptr <= nx
    z(zcnt) = x(ixptr);
    zcnt    = zcnt + 1;
    ixptr   = ixptr + 1;
end
while iyptr <= ny
    z(zcnt) = y(iyptr);
    zcnt    = zcnt + 1;
    iyptr   = iyptr + 1;
end
%--------------------------------------------------------------------
end