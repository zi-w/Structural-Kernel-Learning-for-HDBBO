function [new_curr_grid] = filter_grid(curr_grid, target, threshold_target)

new_curr_grid = [];
threshold = 10e10;

for i = 1:size(curr_grid, 1)
    tmp = threshold_target(curr_grid(i,:));
    if tmp <= threshold
        threshold = tmp;
    end
end
        
for i = 1:size(curr_grid, 1)
    if target(curr_grid(i,:)) <= threshold
        new_curr_grid = [new_curr_grid; curr_grid(i,:)];
    end
end

end