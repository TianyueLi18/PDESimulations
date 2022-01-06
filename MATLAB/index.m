    function index = index(i,j,col)
        if j > col
            error("j > col");
        end
        index = (i-1)*col+j;
    end