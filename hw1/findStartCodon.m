function locations = findStartCodon(inputStr)
    %FINDSTARTCODON Summary of this function goes here
    %   Detailed explanation goes here
    locations = strfind(inputStr, 'ATG');
end

