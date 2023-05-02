%Script: Bonus Homework 1: The Ranked-Choice Vote
%Description: Find the winner of an election with ranked-choice votes.
%Author: Brandon Tran
%UID: 705830462

clc;clearvars; close all;  %Clear vars from previous programs
name = input("Input file of votes to count: ", "s");  %Get input for file
switch(name)  %Switch case based on file name
    case 'votes1.mat'
        votesStruct = load('votes1.mat');   %Load the file and set vars
        votesFile = votesStruct.votes;
        winnerFlag = 0;
        iter = 0;
        [numRows, numCols] = size(votesFile);
        fprintf("                ")  %Print counts and round
        for g=1:numCols
            fprintf('%s    ',pad(int2str(g),2))  
        end
        fprintf("\n")
        while winnerFlag == 0 && iter < numCols   %Loop until winner is found or all votes have been counted
            iter = iter+1;  %Vars to track total count of votes
            totalVotes = zeros(1,numCols);  
            [numRows, numCols] = size(votesFile);
            for i=1:numRows
                if(votesFile(i,1) ~= 0)
                    totalVotes(votesFile(i,1)) = totalVotes(votesFile(i,1))+1;  %Count votes for each candidate                  
                end
                
            end
            fprintf("Round %d Totals: ", iter)
            for g=1:numCols
                fprintf("%s  ", pad(int2str(totalVotes(g)),4))
                if totalVotes(g) == min(totalVotes)
                    index = g;   %Track which candidate has the lowest vote count
                end
            end
            fprintf("\n")
            sortedArr = insertSort(totalVotes);   %Create array of sorted indexes
            winner = totalVotes(sortedArr(length(sortedArr)));
                if max(totalVotes)> numRows/2  %If a candidate has more than half the votes,
                    maximum = totalVotes(1);
                    for i=2:length(totalVotes)
                        if maximum < totalVotes(i)
                            maximum = totalVotes(i);   %Find the candidate
                            index = i;
                        end
                    end
                    winner = index;  %See if there is a winner
                    winnerFlag = 1;
                    break;
                else
                    minimum = 10^12;
                    for i=1:length(totalVotes)
                        if totalVotes(i) < minimum && totalVotes(i) ~= 0
                            minimum = totalVotes(i);   %Find the smallest number of total votes for a candidate not yet eliminated
                            index = i;
                        end
                    end
                    votesFile = removeCandidate(index, votesFile);  %Remove losing candidate
                end
        end
            fprintf("Winning Candidate: "  + winner)
            fprintf("\n===== The Weighted-Choice Vote =====\n")  %Print info for points-based count
            fprintf("             \t ")
            for i = 1:numCols
                fprintf("    %s", pad(int2str(i),2))
            end
            fprintf("\n")
            fprintf("Weighted Points:\t ")
            votesStruct = load(name);  %Reload file
            votesFile = votesStruct.votes;
            points = zeros(1,numCols);  %Vars to track points
            k=1;
            for i=1:numRows
                for j=1:numCols
                    points(votesFile(i,j)) = points(votesFile(i,j)) +(numCols-j+1);  %Calculate points for each candidate
                end
            end
            for n=1:length(points)
                fprintf("%s ", pad(int2str(points(n)),5));  %Print points for each candidate
            end
            maximum = points(1);
            for i=2:length(points)
                if maximum < points(i)
                    maximum = points(i);
                    index = i;
                end
            end
            winner = index;
            fprintf("\n")
            sortedArr = insertSort(points);
            fprintf("Winning Candidate: "  + winner)
        
     case 'votes2.mat'
        votesStruct = load('votes2.mat');  %Load the file and set vars
        votesFile = votesStruct.votes;   
        winnerFlag = 0;
        iter = 0;
        [numRows, numCols] = size(votesFile);
                fprintf("                ")  %Print counts and round
        for g=1:numCols
            fprintf('%s    ',pad(int2str(g),2))  
        end
                fprintf("\n")
        while winnerFlag == 0 && iter < numCols   %Loop until winner is found or all votes have been counted
            iter = iter+1;  %Vars to track total count of votes
            totalVotes = zeros(1,numCols);  
            [numRows, numCols] = size(votesFile);
            for i=1:numRows
                if(votesFile(i,1) ~= 0)
                    totalVotes(votesFile(i,1)) = totalVotes(votesFile(i,1))+1;  %Count votes for each candidate                  
                end
                
            end
            fprintf("Round %d Totals: ", iter)
            for g=1:numCols
                fprintf("%s  ", pad(int2str(totalVotes(g)),4))
                if totalVotes(g) == min(totalVotes)
                    index = g;   %Track which candidate has the lowest vote count
                end
            end
            fprintf("\n")
            sortedArr = insertSort(totalVotes);   %Create array of sorted indexes
            winner = totalVotes(sortedArr(length(sortedArr)));
                if max(totalVotes)> numRows/2  %If a candidate has more than half the votes,
                    maximum = totalVotes(1);
                    for i=2:length(totalVotes)
                        if maximum < totalVotes(i)
                            maximum = totalVotes(i);   %Find the candidate
                            index = i;
                        end
                    end
                    winner = index;  %See if there is a winner
                    winnerFlag = 1;
                    break;
                else
                    minimum = 10^12;
                    for i=1:length(totalVotes)
                        if totalVotes(i) < minimum && totalVotes(i) ~= 0
                            minimum = totalVotes(i);   %Find the smallest number of total votes for a candidate not yet eliminated
                            index = i;
                        end
                    end
                    votesFile = removeCandidate(index, votesFile);  %Remove losing candidate
                end
        end
            fprintf("Winning Candidate: "  + winner)
            fprintf("\n===== The Weighted-Choice Vote =====\n")  %Print info for points-based count
            fprintf("             \t ")
            for i = 1:numCols
                fprintf("    %s", pad(int2str(i),2))
            end
            fprintf("\n")
            fprintf("Weighted Points:\t ")
            votesStruct = load(name);  %Reload file
            votesFile = votesStruct.votes;
            points = zeros(1,numCols);  %Vars to track points
            k=1;
            for i=1:numRows
                for j=1:numCols
                    points(votesFile(i,j)) = points(votesFile(i,j)) +(numCols-j+1);  %Calculate points for each candidate
                end
            end
            for n=1:length(points)
                fprintf("%s ", pad(int2str(points(n)),5));  %Print points for each candidate
            end
            maximum = points(1);
            for i=2:length(points)
                if maximum < points(i)
                    maximum = points(i);
                    index = i;
                end
            end
            winner = index;
            fprintf("\n")
            sortedArr = insertSort(points);
            fprintf("Winning Candidate: "  + winner)
        
        
    otherwise
        if ~exist(name, 'file')
         error(name + ' does not exist');
        else
        votesStruct = load(name);
        votesFile = votesStruct.votes;
    winnerFlag = 0;
        iter = 0;
        [numRows, numCols] = size(votesFile);
                fprintf("                ")  %Print counts and round
        for g=1:numCols
            fprintf('%s    ',pad(int2str(g),2))  
        end
                fprintf("\n")
        while winnerFlag == 0 && iter < numCols   %Loop until winner is found or all votes have been counted
            iter = iter+1;  %Vars to track total count of votes
            totalVotes = zeros(1,numCols);  
            [numRows, numCols] = size(votesFile);
            for i=1:numRows
                if(votesFile(i,1) ~= 0)
                    totalVotes(votesFile(i,1)) = totalVotes(votesFile(i,1))+1;  %Count votes for each candidate                  
                end
                
            end
            fprintf("Round %d Totals: ", iter)
            for g=1:numCols
                fprintf("%s  ", pad(int2str(totalVotes(g)),4))
                if totalVotes(g) == min(totalVotes)
                    index = g;   %Track which candidate has the lowest vote count
                end
            end
            fprintf("\n")
            sortedArr = insertSort(totalVotes);   %Create array of sorted indexes
            winner = totalVotes(sortedArr(length(sortedArr)));
                if max(totalVotes)> numRows/2  %If a candidate has more than half the votes,
                    maximum = totalVotes(1);
                    for i=2:length(totalVotes)
                        if maximum < totalVotes(i)
                            maximum = totalVotes(i);   %Find the candidate
                            index = i;
                        end
                    end
                    winner = index;  %See if there is a winner
                    winnerFlag = 1;
                    break;
                else
                    minimum = 10^12;
                    for i=1:length(totalVotes)
                        if totalVotes(i) < minimum && totalVotes(i) ~= 0
                            minimum = totalVotes(i);   %Find the smallest number of total votes for a candidate not yet eliminated
                            index = i;
                        end
                    end
                    votesFile = removeCandidate(index, votesFile);  %Remove losing candidate
                end
        end
            fprintf("Winning Candidate: "  + winner)
            fprintf("\n===== The Weighted-Choice Vote =====\n")  %Print info for points-based count
            fprintf("             \t ")
            for i = 1:numCols
                fprintf("    %s", pad(int2str(i),2))
            end
            fprintf("\n")
            fprintf("Weighted Points:\t ")
            votesStruct = load(name);  %Reload file
            votesFile = votesStruct.votes;
            points = zeros(1,numCols);  %Vars to track points
            k=1;
            for i=1:numRows
                for j=1:numCols
                    points(votesFile(i,j)) = points(votesFile(i,j)) +(numCols-j+1);  %Calculate points for each candidate
                end
            end
            for n=1:length(points)
                fprintf("%s ", pad(int2str(points(n)),5));  %Print points for each candidate
            end
            maximum = points(1);
            for i=2:length(points)
                if maximum < points(i)
                    maximum = points(i);
                    index = i;
                end
            end
            winner = index;
            fprintf("\n")
            sortedArr = insertSort(points);
            fprintf("Winning Candidate: "  + winner)
        end
end
%====================================================================   
%%  Function - removeCandidate(losingCandidate,votes) 
function votes = removeCandidate(losingCandidate,votes)
    % This function takes the total 2d "votes" matrix, searches for all
    % instances of "losingCandidate", and removes them from the array

    %Error checking for votes 
    [numRows, numCols] = size(votes);  %Set variables
    errorFound = false;
    loserFound = false;
    if(~isnumeric(losingCandidate) || numRows < 1 || numCols < 1 || losingCandidate > numCols) %If parameters are valid
        errorFound = true;
    end
    for i=1:numRows   %For entire votes matrix, check for errors
            for j=1:numCols
                if(votes(i,j) > numCols)  
                    errorFound = true;
                end
                if(losingCandidate == votes(i,j))  %Check if losingCandidate is in votes matrix
                    loserFound =true;
                end
            end
    end
    if(~errorFound && loserFound)  %If no errors found and loser is in matrix
        newVotes = zeros(numRows,numCols);  %Create array of new votes
        for r=1:numRows  %For entire matrix
            k = 1;
            for c=1:numCols
                if votes(r,c) ~= losingCandidate  %If cell is not the losing candidate
                    newVotes(r,k) = votes(r,c);  %Add to new votes matrix
                    k = k+1;
                end
            end
        end
       votes = newVotes;  %Set votes to new Votes
    end 
end

%====================================================================   
%%  Function -  insert(x,index) 
function [x,index] = insert(x,index)
    % x is a sorted colunm vec obtained by applying the insertion process to x.
    % index is a column vec with indices of x
    k=(length(x)-1);
    while x(k+1) < x(k) && k>0  %While current cell in array is larger than previous one
        temp = x(k);  %Swap element
        x(k) = x(k+1);
        x(k+1) = temp;
           
        temp = index(k);  %Swap index
        index(k) = index(k+1);
        index(k+1) = temp;

        k = k-1;  %Decrease iteration by 1
        if k < 2
            break
        end
    end
end

%====================================================================   
%%  Function -  insertSort(x)
function index = insertSort(x)
    %Run insert for each array, each one larger than the previous by 1
   sizeArr = length(x);  
    for k=1:sizeArr-1  %For entire array
        [x(1:k+1), index(1:k+1)] = insert(x(1:k+1), 1:k+1);  %Run insertSort
    end
end