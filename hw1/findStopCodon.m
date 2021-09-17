function locations = findStopCodon(inputStr)
    %FINDSTARTCODON Summary of this function goes here
    %   Detailed explanation goes here
    location1 = strfind(inputStr, 'TAA');
    location2 = strfind(inputStr, 'TAG');
    location3 = strfind(inputStr, 'TGA');
    locations = sort([location1, location2, location3]);
end

