function [U,mesh] = WaterFlowSim(h)
%    h_list = [0.1,0.05,0.025,0.01,0.005];
    x_max = 1;
    x_min = -1;
    y_max = 1;
    y_min = 0;
    num_mesh_x = (x_max - x_min)/h;
    num_mesh_y = (y_max - y_min)/h;
    %Grid is (rows+1)x(col+1)
    %numerical grid: (rows x (col-1);
    %A = [row*(col-1)]
    row = num_mesh_y;
    col = num_mesh_x-1;
    A = zeros(row*col, row*col);
    %TODO: Try sparse matrix
    %finite difference scheme: 2nd order accurate one-sided approximation

    %neumann boundary conditions
    %TODO: verify scheme is second order.
    for j=1:col
        ind = index(1,j,col);
        A(ind,ind) = 3*h./2;
        A(ind,index(2,j,col)) = -2*h;
        A(ind,index(3,j,col)) = h./2;
    end
    
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
    %corners
    ind = index(row,1,col);
    A(ind,index(row-1,1,col)) = 1;
    A(ind,index(row,2,col)) = 1;
    A(ind,ind) = -4;
    ind = index(row,col,col);
    A(ind,index(row-1,col,col)) = 1;
    A(ind,index(row,col-1,col)) = 1;
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
    A = A / (h^2);
    
    %result vector
    fvec = ones(row*col,1)*(-1);
    fvec(1:col) = 0;
    U = A\fvec;
    
    %recreate grid of gutter
    mesh = zeros(num_mesh_y+1, num_mesh_x+1);
    for k=1:row*col
        [i,j] = reverse_index(k,col); %index calculated based on inner box of mesh points, not whole mesh.
        mesh(i,j+1) = U(k);
    end
    
    contour(-1:h:1,1:-h:0,mesh);
    r = num_mesh_y+1;
    c = num_mesh_x+1;
    Z = mesh(:);
    X = zeros(size(Z));
    Y = zeros(size(Z));
    index_cnt = 1;
    for j=1:c %x val
        for i = 1:r %y val
            %in = index(i,j,c);
            X(index_cnt) = -1+h*(j-1);
            Y(index_cnt) = 1-h*(i-1);
            index_cnt=index_cnt+1;
        end
    end
    %scatter3(X,Y,Z);
end