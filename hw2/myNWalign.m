function [score, matrix] = myNWalign(seq1, seq2, match, mismatch, gap)
    % This function performs global alignment between two DNA sequences using
    % the Needleman-Wunsch algorithm (Lecture slide page 55. No trackback
    % needed.)
    % Input:
    %   seq1 and seq2 must be DNA characters (case insensititive)
    %   match score should be positive
    %   mismatch and gap score should be negative (or zero)
    % Returns:
    %   score: score for the best global alignment
    %   matrix: dynamic programming matrix, whose dimension is (1+length(seq1), 1+length(seq2))
    %   note that index in MATLAB starts from 1. As a result, matrix(1, 1)
    %   corresponds to F(0, 0) in lecture slides

    seq1 = upper(seq1);
    seq2 = upper(seq2);

    if (~all(ismember([seq1, seq2], 'ACGT')))
        error('input seq contains non-DNA chars');
    end

    % initialization
    score = 0;
    m = length(seq1);
    n = length(seq2);
    matrix = zeros(m + 1, n + 1);
    seq1 = strcat('\ ', seq1);
    seq2 = strcat('\ ', seq2);
    for i = 1:m + 1
        matrix(i, 1) = (i - 1) * gap;
    end

    for j = 1:n + 1
        matrix(1, j) = (j - 1) * gap;
    end

    % provide your own code to compute matrix and score(no return statement
    % needed in MATLAB).

    for i = 2:m + 1
        for j = 2:n + 1
            if i == 1 & j == 1
                break
            elseif i == 1 & j ~= 1
                left = matrix(i, j - 1);
                up = left;
                up_left = left;
            elseif j == 1
                up = matrix(i - 1, j);
                left = up;
                up_left = up;
            else
                up = matrix(i - 1, j);
                left = matrix(i, j - 1);
                up_left = matrix(i - 1, j - 1);
            end
            if seq1(i) == seq2(j)
                curr = match;
            else
                curr = mismatch;
            end

            opt_list = [up + gap, left + gap, up_left + curr];
            matrix(i, j) = max(opt_list);

            if matrix(i, j) > score
                score = matrix(i, j);
            end

        end

    end
