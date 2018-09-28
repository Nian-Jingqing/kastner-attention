function t = testTensor
data = rand([2 2 3]) + 1;

for ny = 1:size(data, 2)
    data(:, ny, :) = data(:, ny, :) + ny;
end

for nz = 1:size(data, 3)
    data(:, :, nz) = data(:, :, nz) * nz;
end

t = Tensor.Tensor;
t.t = data;
t.t = t.removeMarginalMeans();

data
t.t
