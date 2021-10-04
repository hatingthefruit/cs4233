function [score, matrix] = mySWalign(seq1, seq2, match, mismatch, gap)
% This function performs local alignment between two DNA sequences using
% the Smith-Waterman algorithm. (Lecture slide page 69. No trackback
% needed.)
% Input:
%   seq1 and seq2 must be DNA characters (case insensititive)
%   match score should be positive 
%   mismatch and gap score should be negative 
% Returns: 
%   score: score for the best local alignment
%   matrix: dynamic programming matrix, whose dimension is (1+length(seq1), 1+length(seq2))
%   note that index in MATLAB starts from 1. As a result, matrix(1, 1)
%   corresponds to F(0, 0) in lecture slides 

seq1 = upper(seq1);
seq2 = upper(seq2);
if (~all(ismember([seq1,seq2], 'ACGT')))
    error('input seq contains non-DNA chars');
end

% initialization
score = 0;
m = length(seq1);
n = length(seq2);
matrix = zeros(m+1, n+1);

% complete code below to compute matrix and score (no return statement
% needed in MATLAB).

for i = 1:m
    for j = 1:n
        if i == 1 & j == 1
            break
        elseif i == 1 & j ~= 1
            up = 0;
            left = matrix(i, j-1);
            up_left = 0;
        elseif j == 1
            up = matrix(i-1, j);
            left = 0;
            up_left = 0;
        else
            up = matrix(i-1, j);
            left = matrix(i, j-1);
            up_left = matrix(i-1, j-1);
        end
        if seq1(i) == seq2(j)
            curr = match;
        else 
            curr = mismatch;
        end
        opt_list = [0, up + gap, left + gap, up_left + curr];
        matrix(i,j) = max(opt_list);
        if matrix(i, j) > score
            score = matrix(i, j);
        end
    end
end
