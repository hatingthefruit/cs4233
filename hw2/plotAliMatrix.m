function plotAliMatrix(seq1, seq2, matrix, titleString)
% This function plots the dynamic programming matrix for the (global or
% local) alignment of two sequences.
% Input: 
%   seq1 and seq2 are the two sequences that will be displayed together
%   with the alignment matrix
%   matrix: the alignment matrix that will be displayed. dimension must be
%   (1+length(seq1), 1+length(seq2)
%   titleString: title to be dipslayed above the figure

m = length(seq1);
n = length(seq2);
if (any(size(matrix) ~= [m+1, n+1]))
    error ('matrix size does not match seq length');
end

figure();
imagesc(matrix);
colorbar;
xlabel('seq2');
ylabel('seq1');
set(gca, 'ytick', 2:m+1, 'yticklabels', seq1');
set(gca, 'xtick', 2:n+1, 'xticklabels', seq2');
title(titleString);