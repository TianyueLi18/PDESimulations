function out = WaterFlowSimError
    h_list = [0.1,0.0625,0.05,0.025];
    [~,mesh_exp] = WaterFlowSim(0.0125);
    max_error = zeros(size(h_list));
    for i=1:length(h_list)
        h = h_list(i);
        [~,mesh_emp] = WaterFlowSim(h);
        col_emp = 2/h+1;
        row_emp = 1/h+1;
        jumps = h / 0.0125;
        error = zeros(row_emp*col_emp,1);
        for k=1:row_emp
            for j=1:col_emp
                exp = mesh_exp((k-1)*jumps+1,(j-1)*jumps+1);
                emp = mesh_emp(k,j);
                error(index(k,j,col_emp)) = abs(exp - emp);
            end
        end
        max_error(i) = max(error);
    end
    plot(log(h_list),log(max_error));
end