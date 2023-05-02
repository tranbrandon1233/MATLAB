%Script: Homework 4: The Split-and-Average Problem and Runge-Kutta Methods Applied to Free Fall Problems
%Description: This program finds the average of the two neighboring cells in an array and places it between those elements and finds the weighted average of the three neighbors for problem 1.
%   Additionally, it uses the Runge-Kutta method to ______ for problem 2
%Author: Brandon Tran
%UID: 705830462

clc;clearvars; close all;

problem = input('Enter the problem number (1 or 2): '); %Input problem to display
switch problem  %Switch statement to select problem
    case 1  %Problem 1
       fprintf('Problem 1:\n');
        x = [0,0,1,1];  %Initialize variables
        y=[0,1,0,1];
        w=[-9,-3,3];
        lengthX = length(x);
        err = .001;
        error = 999;
        counter = 0;
        maxIt = 9999;
        figure  %Plot initial array
        plot(x,y);
        title("X vs. Y Before Splitting and Averaging")
        xlabel("X")
        ylabel("Y")
        while (error > .01) &&  (counter < maxIt) %Loop until error is negligible or max number of loops is reached
            xs = splitPoints(x);  %Run functions for x and y
            xa = averagePoints(xs, w);
            ys=splitPoints(y);
            ya=averagePoints(ys,w);
            x=xa;  %Set x and y to output of averagePoints()
            y=ya;
            dx = xa-xs;  %Find difference between output of averagePoints() and output of splitPoints() for x and y
            dy=ya-ys;
            error = abs(max(sqrt(dx.^2+dy.^2)));  %Calculate error
            counter = counter +1;  %Add 1 to counter            
        end
        figure  %plot figure after changes to x and y
        plot(x,y)
        title("X vs. Y After Splitting and Averaging")
        xlabel("X")
        ylabel("Y")
    case 2  %Problem 2
        maxIter=8;  %Set variables
        dt=zeros(maxIter,1);
        dt(1)=1;
        v=zeros(maxIter,1);
        v(1)=0;
       m = 0.5;
       A = 0.276;
       g = 32.174;
       rho = 0.077;
       C_D = 0.47;
       vInf = sqrt((2*m*g)/(rho*C_D*A));
       timeArray = 0:0.001:7;
       for i=1:length(timeArray)   %Iterate for each time step 
           v(i) = vInf*tanh(g*timeArray(i)/vInf);  %Calculate actual velocity for each time step
       end
        err1 = zeros(maxIter*3,1);  %Initialize arrays for errors of each Runge-Kutta method order
        err2 = zeros(maxIter*3,1);
        err3 = zeros(maxIter*3,1);
        err4 = zeros(maxIter*3,1);
       for j = 1:maxIter  %Iterate for all time from 1 to the maximum
           t=0:dt(j):7;  %Initialize time for each time step in dt array
           nt=length(t);   %Set RK and vExact arrays based on t
           RK1=zeros(nt,1);
           RK2=zeros(nt,1);
           RK3=zeros(nt,1);
           RK4=zeros(nt,1);
           vExact=zeros(nt,1);
            timeArray = 0:0.001:7;  %Set time array for vExact
           count = 0;
           for k=2:length(t) %Iterate from 2 for entire length of t
               RK1(k) = RungeKutta(1,g, rho, C_D, A, m, RK1(k-1), dt(j));  %Run RungeKutta() method and assign output to kth value of each RK array
               RK2(k) = RungeKutta(2,g, rho, C_D, A, m, RK2(k-1), dt(j));
               RK3(k) = RungeKutta(3,g, rho, C_D, A, m, RK3(k-1), dt(j));
               RK4(k) = RungeKutta(4,g, rho, C_D, A, m, RK4(k-1), dt(j));
               
               vExact(k)= vInf*tanh(g*timeArray(k)/vInf); %Calculate the actual value of the velocity for time step k

               count = count+1;
               err1((j-1)*3+count) = abs(v(k)-RK1(k));  %Caclulate the errors for each order of the RK method and assign it to err arrays
               err2((j-1)*3+count) = abs(v(k)-RK2(k));
               err3((j-1)*3+count) = abs(v(k)-RK3(k));
               err4((j-1)*3+count) = abs(v(k)-RK4(k));
           end
           if j < 4  %Before the 4th iteration of j
               figure  %Plot all of the RK methods and the actual velocity
               plot(timeArray, v, 0:dt(j):7, RK1, 0:dt(j):7, RK2,0:dt(j):7, RK3,0:dt(j):7, RK4)
               title("Comparison of Actual Velocity And Runge-Kutta Methods")
               xlabel("Time");
               ylabel("Velocity");
           end
            dt(j+1) = dt(j)/2;  %Assign the next value of dt 

       end
% print the results using output function
output(1, err1, dt, maxIter);
output(2, err2, dt, maxIter);
output(3, err3, dt, maxIter);
output(4, err4, dt, maxIter);

end
% 1. 
% splitPoints is a function takes in a vector x and
% returns vector xs, where every odd point in xs is x and
% every even point is the midpoint between the point before and after
% the last point is the mid point between last and first point of x
function [xs] = splitPoints(x)
    % intialize vector for xs
    xs = zeros(2 * length(x), 1);
    for i = 1:length(x)
        % populate the new point with previous ones in every other cell
        xs(2*i - 1) = x(i);        
        % populate the remaining cells at midpoint
        if i ~= length(x)
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
        w = w / sw; %Normalize w
    end
    % initialized xa
    xa = zeros(length(xs)/2,1);
    % iterate through xs
    for i = 1:length(xs)
        % cycle i-1 and i+1 through the end points
        % first element
        if i == 1
            im1 = length(xs);
            im3=2;
        elseif i == length(xs)  %last element
            im3 = 1;
            im1=length(xs)-1;
        else  %all other elements
            im1 = i - 1;
            im3 = i+1;
        end
        % last element
        % figure out ip1
        % weighted 3-point-average
        xa(i) = w(1) * xs(im1) + w(2)*xs(i)+w(3)*xs(im3);
    end
end
%Use the Runge-Kutta methods of the order given by method to estimate vkp1
function [vkp1] = RungeKutta(method,g,rho,C_D,A,m,vk,dt)
    if dt <= 0 || isnumeric(dt) ~= 1  %Check if input is valid
        error("The time step must be a positive number.");
    end
    
    switch(method)  %Select order of RK method based on method parameter
        case 1  %Order 1
                c1=dt*func(vk);
                vkp1 = vk+c1;
        case 2  %Order 2
                c1=dt*func(vk);
                c2=dt*func(vk+0.5*c1);
                vkp1 = vk+c2;
        case 3 %Order 3
            c1=dt*func(vk);
            c2=dt*func(vk+0.5*c1);
            c3=dt*func(vk-c1+2*c2);
            vkp1 = vk+c1/6+2*c2/3+c3/6;
        case 4 %Order 4
            c1=dt*func(vk);
            c2=dt*func(vk+c1/2);
            c3=dt*func(vk+c2/2);
            c4=dt*func(vk+c3);
            vkp1 = vk+c1/6+c2/3+c3/3+c4/6;
        otherwise  %Error in all other cases
            error("Invalid method number.");
    end
    function [out] = func(vk)  %Calclulate function for derivative of velocity
        out = (g-(rho*C_D*A)/(2*m)*vk^2);
    end
end
%Print the error ratios of the Runge-Kutta equations
function output(method, error, dt, maxIt)
    fprintf('                         RK%d\n', method)
    fprintf('                   2s     4s     6s\n')
    
    for j = 2:maxIt %For all elements besides first one after the first one
        % display dt with 6 decimal digits
        fprintf("dt = %.6f: ",dt(j));
        for k = 1: 3
            % print the ratio error((j-1)*3 + k) / error(j*3 + k)
            fprintf("%.4f ", error((j-1)*3 + k) / error(j*3 + k));
        end
        % display new line char
        fprintf('\n');
    end
end