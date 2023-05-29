% given a vector of y values that may have duplicates, return a set of x
% values (about some specified center) that are spread out. This enables
% viewing individual points more clearly (e.g. for overlaying on box plots)
%
% y_values - vector of values to be spread out wherever duplicates are
% present
%
% center_x - center x coordinate about which to spread out the data points
%
% max_offset - maximum absolute offset (e.g., pick 0.5 if you want all data
% points to be between [center_x-0.5,center_x+0.5])
% 
function spread_out_x_values = get_spread_out_x_values(y_values,center_x,max_offset)

spread_out_x_values = center_x*ones(length(y_values),1);

for val_idx = 1:length(y_values)
    % find the number of dumplicate values
    val_i = y_values(val_idx);
    copy_indices = find(y_values==val_i);
    n_copies = length(copy_indices);
    
    % spread the values out based on that duplicate, assuming an
    % imaginary point at the extremes of the offset
    x_offsets = linspace(-max_offset,max_offset,n_copies+2);
    x_offsets([1 end]) = [];
    spread_out_x_values(copy_indices) = center_x + x_offsets;
end