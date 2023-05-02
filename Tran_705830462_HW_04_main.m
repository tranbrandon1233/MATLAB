%Script: Homework 4: The Simple Gravity Pendulum Problem and DNA Analysis
%Description: This program finds the average of the two neighboring cells in an array and places it between those elements and finds the weighted average of the three neighbors for problem 1.
%   Additionally, it finds the number of codons between the start and stop codons of a DNA strand for problem 2
%Author: Brandon Tran
%UID: 705830462

clc;clearvars; close all;
problem = input('Enter the problem number (1 or 2): '); %Input problem to display
switch problem  %Switch statement to select problem
    case 1  %Problem 1
       fprintf('Problem 1:\n');
        x = [0,0,1,1];
        y=[0,1,0,1];
        w=[1/6,1/3,1/2];
        err = 10^-3;
        error = 999;
        counter = 0;
        maxIt = 9999;
        figure
        plot(x,y);
        title("X vs. Y Before Splitting and Averaging")
        xlabel("X")
        ylabel("Y")
        while(error > err ||  counter > maxIt)
            xs = splitPoints(x);
            xa = averagePoints(xs, w);
            ys=splitPoints(y);
            ya=averagePoints(ys,w);
            x=xa;
            y=ya;
            dx = xa-xs;
            dy=ya-ys;
            error = max(sqrt(dx^2+dy^2));
            counter = counter +1;
        end
        figure
        plot(x,y)
        title("X vs. Y Before Splitting and Averaging")
        xlabel("X")
        ylabel("Y")
    case 2  %Problem 2
       
end
% 1. 
% splitPoints is a function takes in a vector x and
% returns vector xs, where every odd point in xs is x and
% every even point is the midpoint between the point before and after
% the last point is the mid point between last and first point of x
function [xs] = splitPoints(x)
    % intialize vector for xs
    xs = zeros(2 * length(x), 1);
    for i = 1:length(x)/2
        % populate the new point with previous ones in every other cell
        xs(2*i - 1) = x(i);        
        % populate the remaining cells at midpoint
        if i ~= length(xs)
            xs(2*i) = (x(i) + x(i+1))/2;
            % the value of xs(2*i) is the average between x(i) and x(i+1)
        else
             xs(2 * i) = (x(i) + x(1))/2;
             %the cell past the last one cycles to the first again, so the values of xs(2 * i) is the average between x(i) and x(1)
        end
    end   
end

% 2.
%This function finds the weighted average of the three neighbors.
function [xa] = averagePoints(xs,w)
    % normalize w
    sw = sum(w);
    % error check for sum of weight and nornalize weight scale
    if sw == 0
        error("Sum of weights cannot be zero.");
    else
        w = w / sw;
    end
    % initialized xa
    % iterate through xs
    for i = 1:length(xs)
        % cycle i-1 and i+1 through the end points
        % first element
        if i == 1
            im1 = length(xs);
            im3=2;
        elseif i == length(xs)
            im3 = 1;
            im1=length(xs)-1;
        else
            im1 = i - 1;
            im3 = i+1;
        end
        % last element
        % figure out ip1
        % weighted 3-point-average
        xa(i) = w(1) * xs(im1) + w(2)*xs(i)+w(3)*xs(im3);
    end
end
