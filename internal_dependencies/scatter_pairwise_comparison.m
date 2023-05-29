%{

Make a scatter plot that enables pairwise comparison of two sets of data
points and also plots ratio lines (i.e., a unit diagonal, a 1.5x line on
both sides, a 2x line on both sides)

group_data (required) - cell array with each cell containing a m-by-2 matrix; each
group will matrix will be plotted as a single color, with the first column
on the x axis and the second column on the y axis

color_data (required) - cell array characters or vectors that indicate the
color to plot each group

xlabel_str (required) - label for x axis

ylabel_str (required) - label for y axis

ZERO_FOR_PLOT (optional) - lower bound for axes

scale_factor_values (optional) - specify the factors to plot reference
lines for; by default, set to 1.5 and 2

%}
function [min_value,max_value] = scatter_pairwise_comparison(group_data,color_data,xlabel_str,ylabel_str,ZERO_FOR_PLOT,group_names,scale_factor_values)

assert(iscell(group_data))

if (~exist('ZERO_FOR_PLOT','var') || isempty(ZERO_FOR_PLOT))
    ZERO_FOR_PLOT = 0.1;
end

if (~exist('group_names','var') || isempty(group_names))
    group_names = arrayfun(@(x) num2str(x,'group %d'), 1:length(group_data), 'UniformOutput', false);
end

if (~exist('scale_factor_values','var') || isempty(scale_factor_values))
    scale_factor_values = [1.5 2];
end

% consolidate x & y labels
label_str = {xlabel_str,ylabel_str};


% determine max value and min value of data
max_value = -Inf;
min_value = Inf;
for group_idx = 1:length(group_data)
    max_value_i = 1.1*max(max(group_data{group_idx}));
    max_value = max(max_value,max_value_i);
    
    min_value_i = min(min(group_data{group_idx}));
    min_value = min(min_value,min_value_i);
end

% Plot the data
for group_idx = 1:length(group_data)
    scatter(group_data{group_idx}(:,1),group_data{group_idx}(:,2),40,color_data{group_idx},'filled'); hold on
end

% Plot and label diagonal line, 1.5x line, and 2x line
ref_lines = {};
ref_lines{end+1,1} = plot([ZERO_FOR_PLOT; max_value],[ZERO_FOR_PLOT; max_value],'-','LineWidth',1,'color','k');
hold on;
for val_idx = 1:length(scale_factor_values)
    scale_factor = scale_factor_values(val_idx);
    ref_lines{end+1,1} = plot([ZERO_FOR_PLOT; max_value],[scale_factor*ZERO_FOR_PLOT; scale_factor*max_value],'-','LineWidth',1/scale_factor,'color',[0.8,0.8,0.8]*(1/scale_factor));
    ref_lines{end+1,1} = plot([scale_factor*ZERO_FOR_PLOT; scale_factor*max_value],[ZERO_FOR_PLOT; max_value],'-','LineWidth',1/scale_factor,'color',[0.8,0.8,0.8]*(1/scale_factor));
    
    PRECISION = 6;
    if (round(mod(scale_factor,1),PRECISION)==0)
        number_format = '%d';
    else
        number_format = '%0.1f';
    end
    text((max_value/scale_factor),(max_value),num2str(scale_factor,number_format),'FontName','Arial','FontSize',14,'HorizontalAlignment','center');
    text((max_value),(max_value/scale_factor),num2str(scale_factor,['1/',number_format]),'FontName','Arial','FontSize',14,'HorizontalAlignment','center');
end

if (max_value~=-Inf)
    xlim([0; max_value]); ylim([0; max_value])
end

% formatting and labeling
set(gca,'FontName','Arial','FontSize',14)
xlabel(label_str{1}); ylabel(label_str{2});
box off
legend(group_names,'Location','Best')

% Make the x axis (and z axis...) length equal to the y axis length
current_axis = gca;
y_axis_length = current_axis.PlotBoxAspectRatio(2);
current_axis.PlotBoxAspectRatio = repmat(y_axis_length,1,3); 

% Move the ref lines to the bottom of the UI stack
cellfun(@(x) uistack(x,'bottom'), ref_lines)

% Print within-group and all-group comparison of differences (median [min, max])
%%% Within-group: For each group, print the within-group differences

% ***** Summarize the raw values for each group
fprintf('***** ***** Summarize the raw values for each group\n');
data_all_groups = [];
for i = 1:length(group_data)
    data_i = group_data{i};
    
    fprintf('\tgroup %d summary (%s): \n',i,group_names{i});
    for j = 1:size(data_i,2)
        print_format = '%0.3f';
        fprintf('\t\tsummary of %s: \n',label_str{j});
        print_summary_stats(data_i(:,j),print_format);
    end
    
    data_all_groups = vertcat(data_all_groups,data_i);
end
fprintf('\tall-group summary:\n');
for j = 1:size(data_all_groups,2)
    print_format = '%0.3f';
    fprintf('\t\tsummary of %s: \n',label_str{j});
    print_summary_stats(data_all_groups(:,j),print_format);
end

% ***** Summarize the ratios of x and y 
fprintf('***** Summarize the ratios (x: %s ; y: %s) \n',label_str{1},label_str{2});
log_data_ratios_all_groups = [];
for i = 1:length(group_data)
    data_i = group_data{i};
    
    data_ratios = data_i(:,2)./data_i(:,1);
    % handle negative ratios
    if any(data_ratios<0)
        warning('taking absolute value for the upcoming ratios since some ratios were negative');
        data_ratios = abs(data_ratios);
    end
    % take log and then later pass exp; this makes inverting ratios for
    % printouts straightforward
    log_data_ratios = log(data_ratios);
        
    fprintf('\tgroup %d summary (%s): \n',i,group_names{i});
    print_format = '%0.3f';
    fprintf('\t\tsummary of y / x: \n');
    print_summary_stats(exp(log_data_ratios),print_format);
    fprintf('\t\tsummary of x / y: \n');
    print_summary_stats(exp(-log_data_ratios),print_format);
        
    log_data_ratios_all_groups = vertcat(log_data_ratios_all_groups,log_data_ratios);
%     fprintf('\t\t%0.6f [%0.6f, %0.6f]\n',median(ratios_y_to_x),min(ratios_y_to_x),max(ratios_y_to_x));
end

%%% All-group: Combine all groups into one; print the whole-group
%%% differences
fprintf('\tall-group differences y / x:\n');
print_summary_stats(exp(log_data_ratios_all_groups),print_format);
fprintf('\tall-group differences x / y:\n');
print_summary_stats(exp(-log_data_ratios_all_groups),print_format);

% print the various summary stats with a given print format
function print_summary_stats(data_ij,print_format)

flag_print_order = 2;

switch flag_print_order
    case 1 % median [min,quant25,quant75,max]
        fprintf('\t\t\tmedian [min,quant25,quant75,max]: %s [%s, %s, %s, %s]\n',...
            num2str(median(data_ij),print_format),...
            num2str(min(data_ij),print_format),...
            num2str(quantile(data_ij,0.25),print_format),...
            num2str(quantile(data_ij,0.75),print_format),...
            num2str(max(data_ij),print_format));
    case 2 % median (quant25, quant75) [min,max]
        fprintf('\t\t\tmedian (quant25, quant75) [min,max]: %s (%s, %s) [%s, %s]\n',...
            num2str(median(data_ij),print_format),...
            num2str(quantile(data_ij,0.25),print_format),...
            num2str(quantile(data_ij,0.75),print_format),...
            num2str(min(data_ij),print_format),...
            num2str(max(data_ij),print_format));
        
    case 3 % median (quant75-quant25 (i.e., IQR)) (quant25, quant75) [min,max]
        fprintf('\t\t\tmedian (IQR) (quant25, quant75) [min,max]: %s (%s) (%s, %s) [%s, %s]\n',...
            num2str(median(data_ij),print_format),...
            num2str(quantile(data_ij,0.75)-quantile(data_ij,0.25),print_format),...
            num2str(quantile(data_ij,0.25),print_format),...
            num2str(quantile(data_ij,0.75),print_format),...
            num2str(min(data_ij),print_format),...
            num2str(max(data_ij),print_format));
end
