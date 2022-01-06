function [i,j] = reverse_index(ind,col)
    if ind < 0
        error("index<0");
    end
    j = mod(ind,col);
    if j == 0
        j = j+ col;
    end
    i = (ind-j) / col+1;
end