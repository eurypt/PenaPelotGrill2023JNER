function fig_handle = compare_group_fitness(fitness_vals_cell_array,fig_handle,ROUND_PRECISION_FOR_SPREAD,colors_cell_array,scatter_size)
if (nargin<2)
    fig_handle = figure('color',[1 1 1 ],'position',[669,98,560,839]);
end

if (nargin<3)
    ROUND_PRECISION_FOR_SPREAD = 2; % for determining how much to spread the values, the numbers are rounded to 10^ROUND_PRECISION_FOR_SPREAD
end

if (nargin<4)
    colors_cell_array = repmat({'k'},length(fitness_vals_cell_array));
end

if (~exist('scatter_size','var') || isempty(scatter_size))
    scatter_size = 15;
end

for i = 1:length(fitness_vals_cell_array)
    bar_width = 0.85;
    boxplot(fitness_vals_cell_array{i},'Positions',i,'widths',bar_width,'colors','k');
    hold on;
    center_x = i;
    max_offset = bar_width/2;
    spread_out_x_values = get_spread_out_x_values(...
        round(fitness_vals_cell_array{i},ROUND_PRECISION_FOR_SPREAD),...
        center_x,max_offset);
    scatter(spread_out_x_values,fitness_vals_cell_array{i},scatter_size,colors_cell_array{i},'filled','LineWidth',1,'MarkerEdgeColor','none')
%     plot(spread_out_x_values,fitness_vals_cell_array{i},'.','MarkerSize',10,'color',colors_cell_array{i})
%     text(spread_out_x_values,fitness_vals_cell_array{i},arrayfun(@(x) num2str(x),1:length(fitness_vals_cell_array{i}),'UniformOutput',false));
    set(gca,'FontSize',12)
    ylabel('rms amplitude (mA)')
end
xlim([0 length(fitness_vals_cell_array)+1]);
min_fitness = min(cellfun(@min, fitness_vals_cell_array));
max_fitness = max(cellfun(@max, fitness_vals_cell_array));
range_fitness = max_fitness-min_fitness;
ylim([min_fitness-0.1*range_fitness,max_fitness+0.1*range_fitness])
set(gca,'XTick',1:length(fitness_vals_cell_array),'XTickLabel',1:length(fitness_vals_cell_array),...
    'XTickLabelRotation',90)

% Make table for KW test
data_array = [];
label_array = [];
for i = 1:length(fitness_vals_cell_array)
    data_array = vertcat(data_array,fitness_vals_cell_array{i});
    label_array = vertcat(label_array,i*ones(length(fitness_vals_cell_array{i}),1));
end

% Run KW test
[p_value,tbl,test_struct] = kruskalwallis(data_array,label_array)


% Swith the active figure
figure(fig_handle)