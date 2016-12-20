function plotPaths(f, indices)

tree = load(f,'-ascii');

[M,J] = size(tree);

if nargin < 2
    indices = (1:M).';
end

figure();
clf;

plot(tree(indices,:));
title('Softbit evolution');
ylabel('Probability of one');
xlabel('Iteration');
end