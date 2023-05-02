%Script: Homework 7: The Game of Life and Bending of an Euler-Bernoulli
%Description: Problem 1 simulates The Game of Life by John Conway, where a cell is killed if it does not have two 
%   or three neighbors or respawns when it has three neighbors, and problem 2 solves a linear system of equations for the 
%   nodal deflections of a simply-supported aluminum beam
%Author: Brandon Tran
%UID: 705830462

clc;clearvars; close all;  %Clear vars from previous programs
rng('shuffle')  %Shuffle rng

problem = input('Enter the problem number (1 or 2): '); %Input problem to display
switch problem  %Switch statement to select problem
    case 1  %Problem 1
       fprintf('Problem 1:\n');
       grid = rand(100,200);  %Set arrays and vars
       nTrials = 500;
       numAlive = 0;
       initialAlive = 0;
       totalAlive = zeros(nTrials+1,1);
       for i=1:100   %Randomly create a living cell in ~1 in every 10 cells  
           for j=1:200
                if grid(i,j) < 0.9
                    grid(i,j) = 0;  
                else
                    grid(i,j) = 1;
                    initialAlive = initialAlive +1;
                end
           end
       end
       totalAlive(1) = initialAlive; %Set the initial value of the total number of living cells
       figure  %Plot grid
       imagesc(grid);
       drawnow
        newGrid = grid;
       for n=1:nTrials  %Iterate for the entire grid for all trials
           for i=1:100 
               for j=1:200
                    if i == 1 %If on the first row
                        im1 = 100;   %Set top row of neighbors to be cells from the bottom row
                    else   %In all other cases,
                        im1 = i-1;   %Set top row of neighbors to be cells from the row below the current one
                    end
                    if i == 100 %If on the last row
                        ip1 = 1;   %Set top row of neighbors to be cells from the bottom row
                    else   %In all other cases,
                        ip1 = i+1;   %Set top row of neighbors to be cells from the row below the current one
                    end
                    if j == 1 %If on the first column
                        jm1 = 200;   %Set left column of neighbors to be cells from the last column
                    else   %In all other cases,
                        jm1 = j-1;   %Set left column of neighbors to be cells from the column to the left of the current one
                    end
                    if j == 200 %If on the last column
                        jp1 = 1;   %Set right column of neighbors to be cells from the first column
                    else   %In all other cases,
                        jp1 = j+1;   %Set right column of neighbors to be cells from the column to the right of the current one
                    end
                    %Calculate number of living neighbors
                    liveNeighbors = newGrid(im1, jm1) + newGrid(im1, j) + newGrid(im1, jp1) + newGrid(ip1, j) + newGrid(ip1, jm1) + newGrid(ip1, jp1) + newGrid(i, jp1) + newGrid(ip1, jm1);
                    if newGrid(i,j)==1  %If cell is alive,
                        if liveNeighbors ~= 2 && liveNeighbors ~= 3 %Kill it if it does not have 2 or 3 neighbors
                            newGrid(i,j) = 0;
                        end
                    else
                        if liveNeighbors == 3  %If it has three neighbors and is not alive, create a new living cell
                            newGrid(i,j) = 1;
                        else
                            newGrid(i,j) = 0;  %Kill it in all other cases
                        end
                    end
               end
           end
       for i = 1:100  
           for j = 1:200
               if(newGrid(i,j) == 1)
                   numAlive = numAlive +1;   %Count number of living cells at end of trial
               end
           end
       end
       totalAlive(n+1) = numAlive; %Assign number of living cells to array
       numAlive = 0;  %Reset number of living cells
       imagesc(newGrid)  %Draw grid
       drawnow
       end
       figure
       plot(1:nTrials+1, totalAlive)  %Plot number of living cells for each trial
       xlabel("Time");
       ylabel("Number of living cells")
       title("Number of Living Cells Over Time")
    case 2  %Problem 2
       n = 50;  %Set variables
       E = 70*10^9;
       L = 2;
       w = 0.05;
       h = 0.07;
        p = 3000;
        d=1.21;
        dx = L/(n-1);
        x = 0:dx:L;
        I = w*h^3/12;
       A = zeros(n-2,n-2);
        for i = 1:n-2   
           for j = 1:n-2
               if i == j
                   A(i,j) = -2;   %Set diagonal in A to -2
               elseif i+1 == j || i-1==j
                    A(i,j) = 1;  %Set two subdiagonals in A to 1
               end   
           end
       end
       M = zeros(n-2,1); 
       index = 1;
       for i = 2:n-1
           M(index) = p*L/2*x(i)-p/2*x(i)^2;  %Fill array with moments for x values from 2 to n-1
           index = index+1;
       end
       b = dx^2*M/(E*I);   %Calculate b vector and y values
        y = A\b;

        deflection = zeros(1,n);   %Create deflection and moment arrays and fill them 
        deflection(2:size(y)+1) = y;
        moment = zeros(1,n);
        moment(2:size(M)+1) = M;
        figure
        subplot(2,1,1)  %Plot moment and deflection vs. x
        plot(moment, '-s')
        xlabel("X")
        ylabel("Moment")
        title("Moment vs. X")
        subplot(2,1,2)
        plot(deflection)
        title("Deflection vs. X")
        xlabel("X")
        ylabel("Deflection")
        
        yExact = -p*d/(24*E*I)*(d^3-2*L*d^2+L^3);  %Calculate the exact value of y

        for i = 1:n-1   
            if x(i) <= d && x(i+1) > d  %Find element of x that is right before d
                slope = (y(i) - y(i-1))/(x(i+1)-x(i));  %Find slope at that point
                yd = y(i-1) + slope*(d-x(i)); %Find y value at that point
                break;
            end
        end
        error = abs((yExact-yd)/yExact)*100;  %Find and print error of approximation vs the actual value
        fprintf("Percentage error: %.4f percent.\n", error);
    otherwise
        %Error message if problem number input is not 1 or 2
        error("Invalid input.");
end