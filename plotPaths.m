function plotPaths(f, indices)
% Plot decoder paths
%
% Input:
%       f       Filepath to decoder debug output file
%       indices Optional vector of bit indices to plot. If not given, all
%               bits will be plotted

tree = load(f,'-ascii');

[~,M] = size(tree);

if nargin < 2
    indices = (1:M).';
end

figure();
clf;
hold all;

plot(tree(:,indices(:)));
title('Softbit evolution');
ylabel('Probability of one');
xlabel('Iteration');



indices_changing = indices(any(diff(sign(tree(:,indices(:))), 1)>0,1));
fprintf('The following bits have changed: ');
fprintf('%d ', indices_changing);
fprintf('\n');
end