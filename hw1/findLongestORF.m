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
    %for index = 0:2
    %    mod_starts = cat(1, mod_starts, starts(mod(starts, 3) == index))
    %    mod_stops = cat(1, mod_starts, stops(mod(stops, 3) == index))
    %end
    mod_starts = {starts(mod(starts, 3) == 0);starts(mod(starts, 3) == 1);starts(mod(starts, 3) == 2);}
    mod_stops = {stops(mod(stops, 3) == 0);stops(mod(stops, 3) == 1);stops(mod(stops, 3) == 2);}
    for ix = 1:3
        for st = 1: length(mod_starts{ix})
            sp = find(mod_stops{ix} > mod_starts{ix}(st), 1);
            if mod_stops{ix}(sp) + 2 - mod_starts{ix}(st) > maxOrf
                maxOrf = mod_stops{ix}(sp) + 3 - mod_starts{ix}(st);
                orf = [mod_starts{ix}(st), mod_stops{ix}(sp) + 2];
            end
        end
    end
    % for st = 1:length(starts)
    %     start_mod = mod(starts(st), 3);
    %     for sp = find(stops > starts(1), 1):length(stops)
    %         if mod(stops(sp), 3) == start_mod
    %             break
    %         end
    %     end
    %     if sp > length(stops)
    %         continue
    %     else
    %         if stops(sp) + 2 - starts(st) > maxOrf
    %             maxOrf = stops(sp) + 3 - starts(st);
    %             orf = [starts(st), stops(sp) + 2];
    %         end
    %     end
    % end
end