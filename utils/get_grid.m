function [newfixhyp] = get_grid(z, fixhyp)

all_cat = unique(z);
fixhyp.xgrid = cell(1,length(all_cat));

for i = 1:length(all_cat)
    indices = (z==all_cat(i));
    curr_discrete = fixhyp.discrete(indices,:);
    
    num_coor = size(curr_discrete, 1);
    load(['utils/grids/' num2str(fixhyp.num_split) '_' num2str(num_coor) '.mat']);
    coors = grids;
    %coors = nmultichoosek(1:fixhyp.num_split, num_coor)
    fixhyp.xgrid{i} = zeros(size(coors,1), num_coor);
    
    for coor_idx = 1:size(coors,1)
        curr_idx = sub2ind(size(curr_discrete), 1:num_coor, coors(coor_idx,:));
        fixhyp.xgrid{i}(coor_idx,:) = curr_discrete(curr_idx);
    end
end

newfixhyp = fixhyp;

end