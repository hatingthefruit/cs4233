function orf = findLongestORF(strand)
%findLongestORG - Description
%
% Syntax: orf = findLongestORG(strand)
%
% Long description
    maxOrf = 0;
    orf = [];
    starts = findStartCodon(strand);
    stops = findStopCodon(strand);
    mod_starts = [];
    mod_stops = [];
    mod_starts = {starts(mod(starts, 3) == 0);starts(mod(starts, 3) == 1);starts(mod(starts, 3) == 2);};
    mod_stops = {stops(mod(stops, 3) == 0);stops(mod(stops, 3) == 1);stops(mod(stops, 3) == 2);};
    for ix = 1:3
        for st = 1: length(mod_starts{ix})
            sp = find(mod_stops{ix} > mod_starts{ix}(st), 1);
            if mod_stops{ix}(sp) + 2 - mod_starts{ix}(st) > maxOrf
                maxOrf = mod_stops{ix}(sp) + 3 - mod_starts{ix}(st);
                orf = [mod_starts{ix}(st), mod_stops{ix}(sp) + 2];
            end
        end
    end
end