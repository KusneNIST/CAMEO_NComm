function [set] = exclude(complete_set, set_to_exclude)

set = complete_set;
keep = true(1,length(complete_set));
[r,c] = size(set_to_exclude);
if(r > c)
    set_to_exclude = set_to_exclude';
end
for i=set_to_exclude
    iexc = find(complete_set == i);
    if(~isempty(iexc))
        keep(iexc) = false;
    end
end

set = complete_set(keep);