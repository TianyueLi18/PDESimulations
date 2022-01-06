function nrm = TDSE(h,k,endtime)
    t_steps = endtime/k;
    x_max = 1;
    x_min = -1;
    num_mesh = (x_max - x_min)/h;
    %Grid is (rows+1)x(col+1)
    %numerical grid: (rows x (col-1);
    %A = [row*(col-1)]
    row = num_mesh-1;
    col = num_mesh-1;
    A = zeros(row*col, row*col);
    
    %five point formulas
    %left+right points
    for i=2:row-1
        %left wall
        ind = index(i,1,col);
        A(ind,index(i+1,1,col)) = 1;
        A(ind,index(i-1,1,col)) = 1;
        A(ind,index(i,2,col)) = 1;
        A(ind,ind) = -4;
        %right wall
        ind = index(i,col,col);
        A(ind,index(i+1,col,col)) = 1;
        A(ind,index(i-1,col,col)) = 1;
        A(ind,index(i,col-1,col)) = 1;
        A(ind,ind) = -4;
    end
    %bottom
    for j=2:col-1
        ind = index(row,j,col);
        A(ind,index(row-1,j,col)) = 1;
        A(ind,index(row,j-1,col)) = 1;
        A(ind,index(row,j+1,col)) = 1;
        A(ind,ind) = -4;
    end
    %top
    for j=2:col-1
        ind = index(1,j,col);
        A(ind,index(2,j,col)) = 1;
        A(ind,index(1,j-1,col)) = 1;
        A(ind,index(1,j+1,col)) = 1;
        A(ind,ind) = -4;
    end
    %corners
    ind = index(row,1,col);
    A(ind,index(row-1,1,col)) = 1;
    A(ind,index(row,2,col)) = 1;
    A(ind,ind) = -4;
    
    ind = index(row,col,col);
    A(ind,index(row-1,col,col)) = 1;
    A(ind,index(row,col-1,col)) = 1;
    A(ind,ind) = -4;
    
    ind = index(1,1,col);
    A(ind,index(2,1,col)) = 1;
    A(ind,index(1,2,col)) = 1;
    A(ind,ind) = -4;
    
    ind = index(1,col,col);
    A(ind,index(2,col,col)) = 1;
    A(ind,index(1,col-1,col)) = 1;
    A(ind,ind) = -4;
    
    %inside points
    for i=2:row-1
        for j=2:col-1
            ind = index(i,j,col);
            A(ind,index(i-1,j,col)) = 1;
            A(ind,index(i+1,j,col)) = 1;
            A(ind,index(i,j-1,col)) = 1;
            A(ind,index(i,j+1,col)) = 1;
            A(ind,ind) = -4;  
        end
    end
    A = -A*1i / (h^2);
    
    %TODO: check fvec
    fvec = zeros(row*col,1);
    for cnt=1:row*col
        [i,j] = reverse_index(cnt,col);
        x = x_min + i*h;
        y = x_min + j*h; %since domain is square
        fvec(cnt) = 15 / 16 * (1-x^2)*(1-y^2);
    end
    %U = A\fvec;
    U = fvec;
    time_U(:,1) = U;
    nrm(:,1) = norm(U).^2*h^2;
   
    %Time pregressions
    %backward Euler for first step
    inverse = inv(eye(size(A))-k*A);
    temp = inverse*U;
    time_U(:,2) = temp;
    nrm(:,2) = norm(temp).^2*h^2;
    
    
    %further time steps uses BDF 2nd order
    for t=3:t_steps
        U_prime = 4*time_U(:,t-1)-time_U(:,t-2);
        inverse = inv(3*eye(size(A))-2*k*A);
        temp = inverse*U_prime;
        time_U(:,t) = temp;
        nrm(:,t) = norm(temp).^2*h^2;
    end
    
end