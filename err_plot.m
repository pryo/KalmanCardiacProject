function err_plot(error_vector,col_index)

surf(reshape(error_vector(:,col_index),sqrt(size(error_vector,1)),...
    sqrt(size(error_vector,1))));
end
