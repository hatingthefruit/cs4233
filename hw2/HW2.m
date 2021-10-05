% script skeleton for HW2

% close all existing figures;
close all;

% clear workspace
clear

% clear commmand window
clc;

%% this cell is for testing myNWalign (Q1)
% complete the myNWalign function and run this cell to test

seq1 = 'atctcta';
seq2 = 'tcata';

% calls the myNWalign function that you define, and output the score as well as alignment
% matrix in the command window
[score_myNW, matrix_mynw] = myNWalign(seq1, seq2, 1, -1, -1)

% FYI, plot the alignment matrix returned from above func call
plotAliMatrix(seq1, seq2, matrix_mynw, 'global alignment matrix');

% FYI, uncomment below to see result from nwalign in MATLAB bioinformatics
% toolbox. The matrix plotted by nwalign is a transpose of myNWalign; hence
% I used seq2 vs seq1. Also, it does not show the first row and first column of
% the matrix returned by myNWalgin. In addition, the two figures look different as they
% use different value to color mappings. You can change the colormap
% by clicking "edit" => "colormap" or use the colormap function ("help colormap").

%[score_nw, ali_nw] = nwalign(seq2, seq1, 'showscore', true, 'alphabet', 'nt', 'scoringmatrix', eye(4)*2-1, 'gapopen', 1, 'extendgap', 1)

%% this cell is for testing mySWalign (Q2)
% complete the mySWalign function and run this cell to test

seq1 = 'atctcta';
seq2 = 'tcact';

% calls the mySWalign function that you define, and output the score as well as alignment
% matrix in the command window

[score, matrix] = mySWalign(seq1, seq2, 1, -1, -1)

% FYI, plot the alignment matrix returned from above func call
plotAliMatrix(seq1, seq2, matrix, 'local alignment matrix');

% FYI, uncomment below to see result from nwalign in MATLAB bioinformatics
% toolbox. The matrix plotted by swalign is a transpose of mySWalign; hence
% I used seq2 vs seq1. Also, it does not show the first row and first column of
% the matrix returned by mySWalgin.

% [score_nwalign, alignment] = swalign(seq2, seq1, 'showscore', true, 'alphabet', 'nt', 'scoringmatrix', eye(4)*2-1, 'gapopen', 1, 'extendgap', 1)

%% This cell is for Q3.

seqL = 200;
nRepeat = 1000;
rand_score = zeros(1, nRepeat);

% 3.1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Insert your code here to compute the alignment scores for 1000 pairs of
% random sequences and save in array rand_score
last_seq = randseq(seqL);

for seqN = 1:nRepeat
    curr_seq = randseq(seqL);
    [rand_score(seqN), ~] = mySWalign(last_seq, curr_seq, 1, -1, -100);
    last_seq = curr_seq;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% display statistics
disp('Alignment statistics from 1000 pairs of random sequences')
disp(sprintf('Min score: %.3g', min(rand_score)))
disp(sprintf('Max score: %.3g', max(rand_score)))
disp(sprintf('Mean score: %.3g', mean(rand_score)))
disp(sprintf('Median score: %.3g', median(rand_score)))

% show the histogram
figure();
hist(rand_score, min(rand_score):max(rand_score));
xlabel('Alignment score');
ylabel('Frequency');
title('Alignment statistics');

% 3.2

thresholds = 11:13; % scoreing thresholds to calcluate probabilities

counts_observed = zeros(1, length(thresholds)); % how many alignment scores are >= each threshold?
p_theory = zeros(1, length(thresholds)); % use extreme value distribution to estimate number of seg
p_theory_divided_by_k = zeros(1, length(thresholds)); % p_theory, but without k

lamda = log(3); % lamda is solved for match = 1, mismatch = -1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Insert your code here to compute counts_observed, p_theory, and
% p_theory_divided_by_k

for num = 1:nRepeat

    if rand_score(num) >= 11
        index = rand_score(num) - 11 + 1;
        counts_observed(index:3) = counts_observed(index:3) + 1;
    end
end

% counts: how many values in rand_score are >= threshold?
% p_theory: assuming k = 0.1, what is the expected number of segments to
% have score >= threshold
k_th = 0.1;
p_theory = (k_th * seqL * seqL)*exp(-lamda*thresholds);
% p_theory_divided_by_k: ignoring k (or assuming k = 1), what is the expected number of segments to
% have score >= threshold
p_theory_divided_by_k = (seqL * seqL)*exp(-lamda*thresholds);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_observed = counts_observed / nRepeat;

figure()
hold on
plot(p_observed, p_theory_divided_by_k, '-+', 'markersize', 15)
plot(p_observed, p_theory, ':o', 'markersize', 15)
legend('theory\_divided\_by\_k', 'theory\_k=0.1');
xlabel('P\_observed');
ylabel('P\_theory\_divided\_by\_k');
hold off

k = mean(p_observed ./ p_theory_divided_by_k)

%% 3.3

seq1 = randseq(seqL); % rand seq of length seqL
mutationRate = 0.3; % mutation rate
mask = randperm(seqL, round(seqL * mutationRate)); % choose random positions to mutate
seq2 = seq1; % seq2 is a copy of seq1
seq2(mask) = randseq(length(mask)); %selected locations are mutated to random chars (may be the same char as in seq1)

score = mySWalign(seq1, seq2, 1, -1, -100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replace the following two lines below with the correct formulas to compute the e-value
% and p-value with the k value estimated from 3.2. Formula for p-value is
% on slide page 172. If the p-value turns out to be 0, use e-value as
% p-value instead.
% Do NOT add semincolon at the end of the line, so that the result will be
% displayed in command window.


eval = (k * seqL * seqL) * exp(-lamda*score)

pval = 1 - exp(-eval);

if pval == 0
    pval = eval
else
    pval
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
