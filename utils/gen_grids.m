function gen_grids(n)

[x1] = ndgrid(1:n);
grids = [x1(:)];
save(['utils/grids/' num2str(n) '_1.mat'], 'grids');

[x1 x2] = ndgrid(1:n);
grids = [x1(:) x2(:)];
save(['utils/grids/' num2str(n) '_2.mat'], 'grids');

[x1 x2 x3] = ndgrid(1:n);
grids = [x1(:) x2(:) x3(:)];
save(['utils/grids/' num2str(n) '_3.mat'], 'grids');

[x1 x2 x3 x4] = ndgrid(1:n);
grids = [x1(:) x2(:) x3(:) x4(:)];
save(['utils/grids/' num2str(n) '_4.mat'], 'grids');

[x1 x2 x3 x4 x5] = ndgrid(1:n);
grids = [x1(:) x2(:) x3(:) x4(:) x5(:)];
save(['utils/grids/' num2str(n) '_5.mat'], 'grids');

