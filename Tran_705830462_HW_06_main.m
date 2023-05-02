%Script: Homework 6: Problem 1 uses a Monte Carlo simulation to find how many people are required in a group before 
%   it is more likely than not that any two of their birthdays in the year 2018 occur in the same calendar week, and 
%   problem 2 models the behavior of two random walkers on a 2D grid of tiles with a wall and tracks the distances of 
%   the two walkers until they collide.
%Description: The Shared Birthday Problem and Random Walk Collisions
%Author: Brandon Tran
%UID: 705830462

clc;clearvars; close all;  %Clear vars from previous programs
rng('shuffle')  %Shuffle rng

problem = input('Enter the problem number (1 or 2): '); %Input problem to display
switch problem  %Switch statement to select problem
    case 1  %Problem 1
       fprintf('Problem 1:\n');
      birthdays = zeros(50,1);  %Declare arrays and vars
      nArray = zeros(1000000,1);
      a = 1;
      for k=1:1000000   %Repeat 10^6 times
          for i=1:50  %Iterate 500 times (more than needed)
              birthdays(i) = randi([1,365]);  %Create random birthday and store in birthdays array
              bdayFound = false;
              for n=1:i-1  %Repeat for entire filled length of birthdays so far
                  if ceil(birthdays(i)/7) == ceil(birthdays(n)/7) || (birthdays(i) == 365 && ceil(birthdays(n)/7) == 1) || (birthdays(n) == 365 && ceil(birthdays(i)/7) == 1) %If birthdays are on the same week,
                      nArray(a) = i;  %Store the size of the birthdays array in nArray
                      a = a+1;  %Iterate for next index of nArray
                      bdayFound = true;
                      break;
                  end
              end
              if bdayFound %If two birthdays on the same week were found, break from the loop
                break;
              end
          end
            %Reset birthdays array
            birthdays = zeros(50,1);
      end
      if median(nArray) < 10
        fprintf("Median number of people: 0%d\n", median(nArray));
      else
        fprintf("Median number of people: %d\n", median(nArray));
      end

        histogram(nArray, max(nArray)-min(nArray))
        title("Frequency of Number of Birthdays")
        xlabel("Number of Birthdays")
        ylabel("Frequency")
    case 2  %Problem 2
        n = input('\nPlease input the size of the array: ');
        while isnan(n) || mod(n,1) ~= 0 || mod(n, 2) == 0 || n < 11
            fprintf('\nError: number of steps must be a positive integer, odd, and equal to or greater than 11.\n');
            n = input('\nPlease enter the number of steps: ');
        end
               %   Set initial position of particle A
        xA = -n/2+0.5;
        xAk = xA;
        yA = 0;
        yAk = yA;
        %   Set initial position of particle B
        xB = n/2-0.5;
        xBk = xB;
        yB = -n/2+0.5;
        yBk = yB;
        %   Set boundary conditions [up,down,left,right]
        BD = [n/2-0.5,-n/2+0.5,-n/2+0.5,n/2-0.5];
        % Set value of disturbance
        d = rand*0.1-0.05;
       
        %   Set collision flag and initial step to zero
        collisionFlag = false;
        k = 0;
        % Plot particle movements
        figure(1)
        hold on;
        set(gca,'xtick',-n/2:n/2);
        set(gca,'ytick',-n/2:n/2);
        grid on;
        xlim([-n/2 n/2]);
        ylim([-n/2 n/2]);
        axis square;
        title('2D Random Walk','FontSize',24);
        set(gcf,'Position',[30 350 850 450]);
        set(gca,'LineWidth',1,'FontSize',12);
        set(gcf,'renderer','painters');
        print('-fillpage','dpdf')
        % Create value of the y-axis of the wall
        yWall = randi(2*n)-n;
        % Create arrays for locations of A and B and distance between them
        locationsAX = zeros(500,1);
        locationsAY = zeros(500,1);
        locationsBX = zeros(500,1);
        locationsBY = zeros(500,1);
        distances = zeros(500,1);
        % Change value of yWall until it is not within two tiles of boundary
        while yWall > double(n/2)-2 || yWall < double(-n/2)+2 
            yWall = randi(2*n)-n;
        end
        % Set location of yWall
        yWallX1 = [3.5,.5,0.5,3.5];
        yWallX2 = [-3.5,-.5,-0.5,-3.5];
        yWallY1 = [yWall-0.5,yWall-0.5,yWall+0.5,yWall+0.5];
        % Draw yWall using the color black
        DrawRect(yWallX1,yWallY1,0,0,'k')
        DrawRect(yWallX2,yWallY1,0,0,'k')
        %   Iterate while collision has not occured and maximum step did not exceed
        while collisionFlag == false && k < 500
            pause(0.1);

            % Add the locations of A and B and their distances from each other to their respective arrays
            locationsAX(k+1) = xAk;
            locationsAY(k+1) = yAk;
            locationsBX(k+1) = xBk;
            locationsBY(k+1) = yBk;
            distances(k+1) = sqrt((locationsBX(k+1)-locationsAX(k+1))^2+((locationsBY(k+1)-locationsAY(k+1)))^2);

            % Set next location of x and y
            [xAkp1,yAkp1] = RandWalk_2D(xAk,yAk,BD,d,yWall);
            [xBkp1,yBkp1] = RandWalk_2D(xBk,yBk,BD,d,yWall);
            
            %   Create particle A on grid for step (k)
            xAk_val = [xAk - 0.5, xAk + 0.5, xAk + 0.5, xAk - 0.5];
            yAk_val = [yAk - 0.5, yAk - 0.5, yAk + 0.5, yAk + 0.5];
            
            %   Create particle A on grid for step (k+1)
            xAkp1_val = [xAkp1 - 0.5, xAkp1 + 0.5, xAkp1 + 0.5, xAkp1 - 0.5];
            yAkp1_val = [yAkp1 - 0.5, yAkp1 - 0.5, yAkp1 + 0.5, yAkp1 + 0.5];
            
            %   Create particle B on grid for step (k)
            xBk_val = [xBk - 0.5, xBk + 0.5, xBk + 0.5, xBk - 0.5];
            yBk_val = [yBk - 0.5, yBk - 0.5, yBk + 0.5, yBk + 0.5];
            
            %   Create particle B on grid for step (k+1)
            xBkp1_val = [xBkp1 - 0.5, xBkp1 + 0.5, xBkp1 + 0.5, xBkp1 - 0.5];
            yBkp1_val = [yBkp1 - 0.5, yBkp1 - 0.5, yBkp1 + 0.5, yBkp1 + 0.5];
            
    
            %   Redraw figure with "red" as the step (k) particle A,
            %   "yellow" as the step (k) particle B,
            %   "blue" as the step (k+1) particle A,
            %   "green" as the step (k+1) particle B
            
            % Draw particles and their trails
            DrawRect(xAk_val,yAk_val, 0, 0,'r');
            DrawRect(xAkp1_val,yAkp1_val, 0, 0,'b');
            DrawRect(xBk_val,yBk_val,0,0,'y');
            DrawRect(xBkp1_val,yBkp1_val,0,0,'g');
            
            %   Update the new position for step (k+1)
            xAk = xAkp1;  yAk = yAkp1;
            xBk = xBkp1;  yBk = yBkp1;
            %Increase k by 1
            k = k + 1;
            if xAk == xBk && yAk == yBk
                collisionFlag = true;
            end
        end 
        %Plot the distances between A and B
        figure
        plot(1:k, distances(1:k))
        xlabel("Iterations")
        ylabel("Distance")
        title("Distance vs. Iteration")
    otherwise
        %Error message if problem number input is not 1 or 2
        error("Invalid input.");
end
%====================================================================   
%%  Function - RandWalk_2D(x0,y0,BD,d,yWall)
function [x,y] = RandWalk_2D(x0,y0,BD, d, yWall)
%   2D Random Walk Function that takes in initial positions x0 and y0
%   Computes a directional movement based on a random number with a
%       disturbance d for one of the directions
%   Accounts for boundary conditions vector BD [up,down,left,right]
%   Adds a 7x1 wall named yWall centered at x = 0 and with a random
%       y-coordinate that is at least two tiles away from the top and bottom
%       walls of the grid and a hole at the center
    %   Create a random number for calculating movement
    r = rand;
    %   Move North with disturbance d
    if r <= 0.2-d && (y0+1 ~= yWall || (y0+1 == yWall && (x0 < -3.5 || x0 > 3.5)) || x0==0)
        x = x0;
        y = y0 + 1;
        %   Account for upper limit in boundary
        if y >= BD(1)
            y = BD(1);
        end
    %   Move South
    elseif 0.2 < r && r <= 0.4 && (y0-1 ~= yWall || (y0-1 == yWall && (x0 < -3.5 || x0 > 3.5)) || x0==0)
        x = x0;
        y = y0 - 1;
        %   Account for lower limit in boundary
        if y <= BD(2)
            y = BD(2);
        end
    %   Move West
    elseif 0.4 < r && r <= 0.6 && (y0 ~= yWall || (y0 == yWall && (x0 -1 < -3.5 || x0-1 > 3.5)) || x0-1==0)
        x = x0 - 1;
        y = y0;
        %   Account for left limit in boundary
        if x <= BD(3)
            x = BD(3);
        end
    %   Move E
    elseif 0.6 < r && r <= 0.8 && (y0 ~= yWall || (y0 == yWall && (x0 +1 < -3.5 || x0+1 > 3.5)) || x0+1==0)
        x = x0 + 1;
        y = y0;
        %   Account for right limit in boundary
        if x >= BD(4)
            x = BD(4);
        end
    %   No movement
    else
        x = x0;
        y = y0;
    end
end

%====================================================================   
%%  Function - DrawRect(a,b,L,W,c)
function DrawRect(a,b,L,W,c)
% Function that draws each colored tile as needed on the x-coordinate
% given by a, y-coordinate given by b, displacement of each based on L and
% W respectively, and using the color given by c
    x = [a a+L a+L a];
    y = [b b b+W b+W];
    fill(x,y,c);
end