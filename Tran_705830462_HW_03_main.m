%Script: Homework 3: The Three Species Problem and the Pocket Change
%Description: This program calculates the population of three different
%   species overtime based on given differential equations for problem 1.
%   It also calculates the average number of coins to create change from 1
%   to 99 cents.
%Author: Brandon Tran
%UID: 705830462

clc;clearvars; close all;
problem = input('Enter the problem number (1 or 2): '); %Input problem to display
switch problem  %Switch statement to select problem
    case 1  %Problem 1
       fprintf('Problem 1:\n');
       stepSize = 10^-6;
       x0 = 10;
       y0 = 5;
       z0 = 2.5;
       tol = 10^-12;
       tFinal = 20;
       timeStep = ceil(tFinal/stepSize);
       disp("Time      X       Y       Z       it");
       for it = 0:timeStep
        x = x0 + stepSize*(x0*(1-x0/10)-0.75*x0*y0-2*x0*z0);
        y = y0 + stepSize*(1.5*y0*(1-y0/5)-0.5*y0*x0-1.5*y0*z0);
        z = z0+stepSize*(3*z0*(1-z0/2.5)-1.5*z0*x0-0.5*z0*y0);
        diffX = x - x0;
        diffY = y-y0;
        diffZ = z-z0;
        diffCount = 0;
        while abs(diffX) > tol || abs(diffY) > tol || abs(diffZ) > tol
            newX = x0 + stepSize*(x*(1-x/10)-0.75*x*y-2*x*z);
            newY = y0 + stepSize*(1.5*y*(1-y/5)-0.5*y*x-1.5*y*z);
            newZ = z0+ stepSize*(3*z*(1-z/2.5)-1.5*z*x-0.5*z*y);
            diffX = x - newX;
            diffY = y - newY;
            diffZ = z - newZ;
            x = newX; 
            y = newY; 
            z = newZ; 
            diffCount = diffCount +1;
        end
        if mod(it*stepSize,1) == 0 
            fprintf("   %d  %.2d  %.2d %.2d  %d\n", it*stepSize,x,y,z,diffCount);
        end
        x0 = x;
        y0 = y;
        z0 = z;
       end
    case 2  %Problem 2
        fprintf('Problem 2:\n')
        coinCount = 0;
        for n = 1:99
            totalChange = n;
            coinCount = floor(totalChange/25);
            totalChange = mod(totalChange, 25);
            coinCount = coinCount + floor(totalChange/10);
            totalChange = mod(totalChange, 10);
            coinCount = coinCount + floor(totalChange/5);
            totalChange = mod(totalChange, 5);
            coinCount = coinCount + totalChange;
            fprintf("Average Number of Coins = %.2f\n", coinCount)
        end
end
   
