%Script: Bonus Homework 2: Newton's Method
%Description:  Applies both Newton's method and a modified version of Newton's
%   method to find the roots of an arbitrary function f iteratively in terms of various initial
%   guesses.
%Author: Brandon Tran
%UID: 705830462

clc;clearvars; close all;  %Clear vars from previous programs
rng('shuffle');

method = input("Input the method number: ");  %Prompt user to input method number
eps = 10^-12;  %Set variables
iterMax = 1000;
x0 = -4:0.01:4;
xk = zeros(length(x0), 1);
iter = zeros(length(x0), 1);
f0= zeros(length(x0), 1);
f = @(x) sin(x)*cos(2*x)-x^2/10+1;
derivF = @(x) cos(x)*cos(2*x) - 2*sin(x)*sin(2*x)-x/5;
for i=1:length(x0)  %Iterate for all of x0
    f0(i) = f(x0(i));   %Set f0 to function for value of x0
     [xk(i),iter(i)] = NewtonsMethod(f,derivF,x0(i),eps,iterMax,method);  %Use NewtonsMethod to solve
end
figure  %Plot x0 vs f0 and xk
plot(x0,f0)
hold on
plot(x0,xk)
title("x vs. f0 and x vs. xk")
xlabel("x")
ylabel("f0 and xk")
legend("f0","xk")

hold off

figure  %Plot bar chart for number of iterations
bar(x0, iter)
title("Number of Iterations for Each Value of x0")
xlabel("Value of x0")
ylabel("Number of iterations")

function [xk,iter] = NewtonsMethod(f,derivF,x0,eps,iterMax,method)
%Performs the Newton and modified Newton methods to find the roots of a
%function iteratively in terms of various initial guesses
%Parameters: f is the function, derivF is the derivative of the function,
%x0 is the initial guess, eps is the tolerance, iterMax is the maximum
%number of iterations, and method is the method chosen to find the roots
    if ~isscalar(x0) || ~isscalar(eps) || ~isscalar(iterMax) || ~isscalar(method) || ~isnumeric(x0) || ~isnumeric(eps) || ~isnumeric(iterMax) || ~isreal(eps) || ~isreal(iterMax) || eps <= 0 || iterMax < 1 || rem(iterMax,1) ~= 0
        error("Invalid parameter(s) for NewtonsMethod function") %Error message if an invalid parameter is found
    end
    iter = 0;  %Initialzie variables
    xk = x0;
    h = 0.5;
    while abs(f(xk)) > eps && iter < iterMax  %Loop until iterMax is reached or function is below tolerance
        if method == 1
            xkp1=xk-f(xk)/derivF(xk);  %Method 1
        else
            derivFversion1 = (f(xk+h)-f(xk))/h;  %Both versions of method 2
            derivFversion2 = (f(xk)-f(xk-h))/h;

            xkp1Version1 = xk-f(xk)/derivFversion1;
            xkp1Version2 = xk-f(xk)/derivFversion2;

            if (xkp1Version2 > xkp1Version1)  %Select the smaller one
                xkp1  = xkp1Version1;
            else 
                xkp1 = xkp1Version2;
            end
            h = abs(xkp1-xk)/2;  %Update h
        end
        xk = xkp1;   %Update vars
        iter = iter +1;
    end
end