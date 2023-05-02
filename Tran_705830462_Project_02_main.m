%Script: Final Project
%Description:  The goal of this project is to find the unknown temperature distribution of an
% aluminum fin that is subjected to a heat input using the point
% collocation method and the radial basis function as the interpolant
%Author: Brandon Tran
%UID: 705830462

clc;clearvars; close all;  %Clear vars from previous programs
rng('shuffle');

function [nNodes,nodes] = defineNodes(nNodesX,nNodesY,tol,L,H, whichPart,random)
%Create and sort nodes from 0 to L
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
                    nodes(i) = tol + (L - 2*tol) * rand(1, nNodesX); %Random;y redefine nodes
                    nodes(2:length(nodes)-1) = mergeSort(nodes(2:length(nodes)-1)); %Sort nodes again
                    nodeMoved = false;
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
        end
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