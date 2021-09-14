function revComp = getReverseComp(dna_str)
    keys = {'A', 'C', 'G', 'T'};
    values = {'T','G', 'C', 'A'};
    valMap = containers.Map(keys, values);
    x = 1;
    len = strlength(dna_str);
    revComp = dna_str;
    while x <= len
        revComp(x) = valMap(dna_str(x));
        x = x +1;
    end
    revComp = reverse(revComp);
end
%% modify this code so that revComp is the reverse complementary of dna_str