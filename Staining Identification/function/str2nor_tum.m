%% Divide normal stroma and tumor stroma according to the category of adjacent tissues
function result = str2nor_tum(str, mask_new, side, prop)

if str == 3
    blockproc = mask_new(side(1):side(2), side(3):side(4));
    
    if (isempty(blockproc(blockproc(:) == 2))) && (isempty(blockproc(blockproc(:) == 4)))
        if prop(2) >= prop(4)
            result = 3;
        else
            result = 5;
        end
        
    else
        if length(blockproc(blockproc(:) == 2)) >= length(blockproc(blockproc(:) == 4))
            result = 3;
        else
            result = 5;
        end
    end
    
else
    result = str;
end

end