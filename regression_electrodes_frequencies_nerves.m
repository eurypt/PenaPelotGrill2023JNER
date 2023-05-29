%{

Plot Figures 3, 4, and 5 quantifying the effects of nerve, electrode, and
frequency on block threshold current, block threshold power, and onset
response.

%}

function regression_electrodes_frequencies_nerves()

analyze('frequency_tests_data_table_max.mat','current');
analyze('frequency_tests_data_table_max.mat','power');
analyze('frequency_tests_data_table_max.mat','onset');

% metric_str (required) - set to 'current' or 'power' or 'onset'
function analyze(filename,metric_str)

% add paths to external dependencies (i.e., functions I obtained from
% Mathworks File Exchange for plotting)
addpath('external_dependencies/shadedErrorBar');
addpath('external_dependencies/subaxis');

% add path to internal dependencies (i.e., functions I wrote to facilitate
% plotting/analysis/etc.)
addpath('internal_dependencies');

assert(any(strcmp(metric_str,{'current','power','onset'})),...
    'metric_str must be set to current or power or onset');

load(filename,'frequency_tests_data_table')

% specify how to combine replicates; 'mean' takes the mean; 'median' takes
% the median; 'alt_median' takes the value closest to the median unless two
% values are tied for closest in which case take the smaller of the two
average_type = 'median';

% set to true to use a sum of terms rather than a product of terms ;
regression_parameters.flag_log_transform = true;

% set this to exp(1), 10, 2, or whatever the desired log base transform is
regression_parameters.log_base = 2;

regression_parameters.include_random_effects = false;

regression_parameters.include_electrode_effects = true;

% set to 'categorial' or 'continuous'
regression_parameters.frequency_type =  'categorical'; %'categorical'; %'continuous';

flag_include_hybrid_interactions = 1;

freq_cutoff = 14;

% max number of times to reweight
reweight_count_max = 10;

% set to true to peform weighted least squares regression ; set to false
% otherwise; the weights will be based on the frequency
regression_parameters.weighting = false;

regression_parameters.reweight = false;


%%% Check inputs
assert(strcmp(regression_parameters.frequency_type,'categorical')||strcmp(regression_parameters.frequency_type,'continuous'),...
    'Error: regression_parameters.frequency_type should be set to either ''categorical'' or ''continuous''')

%%% Define log transform func and inverse
log_transform_func = @(x) log(x)./log(regression_parameters.log_base);
inv_log_transform_func = @(x) regression_parameters.log_base.^(x);


%% Section 5: Select metric, truncate expdate, remove replicates

%%% Remove rows that do not have a valid metric value
if (strcmp(metric_str,'onset'))
    frequency_tests_data_table(isnan(frequency_tests_data_table.onset_response_at_single_amp_Ns),:) = [];
end

%%% Select the metric that will be used for subsequent plotting and
%%% analysis
metric = zeros(size(frequency_tests_data_table,1),1);
if (strcmp(metric_str,'current'))
    metric = frequency_tests_data_table.thresholds_amplitude_mA;
    thresholds_str = 'thresh50 (mA)';
elseif (strcmp(metric_str,'power'))
    metric = frequency_tests_data_table.thresholds_amplitude_mW;
    thresholds_str = 'thresh50 (mW)';
elseif (strcmp(metric_str,'onset'))
    metric = frequency_tests_data_table.onset_response_at_single_amp_Ns;
    thresholds_str = 'onset response (N*s)';
end

%%% Identify rows which have replicates, then delete all but one of those
%%% rows, and store the average of the replicates into that row
% find all the unique expdate and frequency combos
all_expdate_freq_combos = [str2double(frequency_tests_data_table.expdate),frequency_tests_data_table.frequencies_kHz];
unique_expdate_freq_combos = unique(all_expdate_freq_combos,'rows');

% for each unique combo...
for i = 1:length(unique_expdate_freq_combos)
    % find all indices of rows having those combo values
    all_indices_of_combo_i = find(ismember(all_expdate_freq_combos,unique_expdate_freq_combos(i,:),'rows'));
    assert(length(all_indices_of_combo_i)>=1,'expected every combo of expdate and freq to have at least one matching row...')
    
    % if are multiple rows having this combo...
    if (length(all_indices_of_combo_i)>1)
        % set the metric value of first row having that combo to the average of
        % the metric values
        if (strcmp(average_type,'mean'))
            
            % Identify the latest test and keep only the row in
            % frequency_tests_data_table corresponding to that test; keep only
            % that row
            [~,datetime_sortorder] = sort(frequency_tests_data_table.datetime_adicht_last_saved(all_indices_of_combo_i));
            [~,max_idx] = max(datetime_sortorder);
            idx_of_row_to_remove = all_indices_of_combo_i(setdiff(1:length(all_indices_of_combo_i),max_idx));
            
            % For metric to plot, take the mean of the replicates
            metric(all_indices_of_combo_i(max_idx)) = ...
                mean(metric(all_indices_of_combo_i));
            
            % Remove the rows identified to remove; remove the same row
            % from both the metric variable and the frequency_tests_data_table
            frequency_tests_data_table(idx_of_row_to_remove,:) = [];
            metric(idx_of_row_to_remove,:) = [];
        elseif (strcmp(average_type,'median'))
            
            % Identify the latest test and keep only the row in
            % frequency_tests_data_table corresponding to that test; keep only
            % that row
            [~,datetime_sortorder] = sort(frequency_tests_data_table.datetime_adicht_last_saved(all_indices_of_combo_i));
            [~,max_idx] = max(datetime_sortorder);
            idx_of_row_to_remove = all_indices_of_combo_i(setdiff(1:length(all_indices_of_combo_i),max_idx));
            
            % For metric to plot, take the median of the replicates
            metric(all_indices_of_combo_i(max_idx)) = ...
                median(metric(all_indices_of_combo_i));
            
            % Remove the rows identified to remove; remove the same row
            % from both the metric variable and the frequency_tests_data_table
            frequency_tests_data_table(idx_of_row_to_remove,:) = [];
            metric(idx_of_row_to_remove,:) = [];
            
        elseif (strcmp(average_type,'alt_median'))
            % keep only the value that is closest to the median of the
            % replicates; delete all others; this behavior reduces to
            % simply taking the median when the number of replicates is odd
            [~,idx_of_closest_to_mean] = min(abs(metric(all_indices_of_combo_i)-median(metric(all_indices_of_combo_i))));
            
            % delete all other rows
            frequency_tests_data_table(all_indices_of_combo_i(setdiff(1:length(all_indices_of_combo_i),idx_of_closest_to_mean)),:) = [];
            metric(all_indices_of_combo_i(setdiff(1:length(all_indices_of_combo_i),idx_of_closest_to_mean))) = [];
        else
            error('Invalid value of average_type; must be either mean, median, or alt_median')
        end
        
        % recalulate the all_expdate_freq_combos matrix because the number
        % of rows in frequency_tests_data_table changed during the deletion
        % above
        all_expdate_freq_combos = [str2double(frequency_tests_data_table.expdate),frequency_tests_data_table.frequencies_kHz];
    end
end

assert(size(metric,1)==size(frequency_tests_data_table,1),'metric and frequency_tests_data_table should have both had the same replicates removed')


%% Section 6: Fit polynomials to threshold vs. frequency


if (~strcmp(metric_str,'onset'))
    % Define inputs
    flag_weight_polynomial_fits = 0; % set to 1 to weight the measurements based on their value; this establishes that uncertainty in block threshold is greater at higher frequencies
    flag_iterate_nerve_effects = 0; %  set to 1 to scale the values of each nerve to bring all the measured values closer to the mean, thus reducing variability in the data while still taking advantage of the repeated measures within nerves
    std_type = 'mean'; % set to 'mean' to plot std of the mean of the fits (as opposed to std of samples ('sample'))
    
    % Then, for each nerve, we identified a scaling factor that minimized the
    % sum of squared error between the prediction and the measured data, and
    % we iteratively recalculated the polynomial fits until the sum of
    % squared error between the predictions and all the data exhibited
    % negligible change.
    
    format compact
    
    cuff_types = {'monopolar','bipolar','tripolar'};
    
    if (flag_weight_polynomial_fits)
        weighting_function_polynomial_fits = @(y_vals) (1./y_vals).^2;
        flag_scale_std_by_y = true;
    else
        weighting_function_polynomial_fits = @(y_vals) ones(size(y_vals));
        flag_scale_std_by_y = false;
    end
    
    for i=1:length(cuff_types)
        indices_i = strcmp(frequency_tests_data_table.electrode_type,cuff_types{i});
        x = frequency_tests_data_table.frequencies_kHz(indices_i);
        x_shift = 20; % [kHz] shift by 20 kHz
        y = metric(strcmp(frequency_tests_data_table.electrode_type,cuff_types{i}));
        
        % initialize scaling factors to 1
        scaling_factors = ones(size(y));
        
        % if flag_iterate_nerve_effects == 1, set a small tolerance to iterate
        % scaling factors; otherwise, set Infinite tolerance to *not* iterate
        % at all
        if (flag_iterate_nerve_effects)
            tol = 0.01; % tolerance; when gof SSE changes by less than this amount, stop iterating the scaling factor calculations
        else
            tol = Inf;
        end
        
        % Perform fitting and (while tolerance has not been reached) iterate
        % across scaling factor values
        delta_gof_sse = Inf;
        latest_gof_sse = Inf;
        iterate = true;
        while (iterate)
            [fit_obj,gof] = fit(x-x_shift,(y.*scaling_factors),'poly2','Weights',weighting_function_polynomial_fits(y.*scaling_factors)); % Same as the following manual process: INPUT_MATRIX = [((x-x_shift).^2)./y, (x-x_shift)./y, 1./y]; OUTPUT_MATRIX = ones(length(y),1); fit_coeffs = (INPUT_MATRIX\OUTPUT_MATRIX);
            
            if (latest_gof_sse~=Inf)
                delta_gof_sse = (gof.sse-latest_gof_sse)/latest_gof_sse;
            end
            latest_gof_sse = gof.sse;
            
            if (delta_gof_sse<=tol)
                iterate = false;
            else
                y_predicted = feval(fit_obj,x-x_shift);
                [unique_expdate_cuff_type_i,~,IC] = unique(frequency_tests_data_table.expdate(indices_i));
                unique_scaling_factors = cellfun(@(input) y(strcmp(frequency_tests_data_table.expdate(indices_i),input))\...
                    y_predicted(strcmp(frequency_tests_data_table.expdate(indices_i),input)),...
                    unique_expdate_cuff_type_i);
                scaling_factors = unique_scaling_factors(IC);
            end
        end
        
        
        fprintf('\ncuff type: %s\n',cuff_types{i})
        flag_plot_scatter_data = 0; % set to 0 to *not* plot raw xy data onto the fit
        if (i==1)
            fig_handle = figure('color',[1 1 1]);
        end
        fit_and_plot_xy_data(x,y.*scaling_factors,x_shift,fit_obj,gof,std_type,flag_scale_std_by_y,flag_plot_scatter_data,fig_handle)
        xlabel('frequency (kHz)');
        ylabel(thresholds_str);
        %     title(cuff_types{i})
        
        % set figure size as desired
        set(gcf,'position',[ 231         552        1085         418])
    end
    
    % change line width of all lines
    gg=get(gca,'Children');
    for i=1:length(cuff_types)
        gg(i).LineWidth=2;
    end
    
    % make legend
    legend(cuff_types,'Location','NorthWest')
    
    % set xlim to go to zero
    xlim_i = get(gca,'XLim');
    xlim([0 xlim_i(2)])
    %{
Incomplete / Archived / Obsolete


    %{
for i=1:length(cuff_types)
    x = frequency_tests_data_table.frequencies_kHz(strcmp(frequency_tests_data_table.electrode_type,cuff_types{i}));
    x_shift = 20; % [kHz] shift by 20 kHz
    y = metric(strcmp(frequency_tests_data_table.electrode_type,cuff_types{i}));
    [fit_obj,gof] = fit(x-x_shift,y,'poly2','Weights',weighting_function_polynomial_fits(y)); % Same as the following manual process: INPUT_MATRIX = [((x-x_shift).^2)./y, (x-x_shift)./y, 1./y]; OUTPUT_MATRIX = ones(length(y),1); fit_coeffs = (INPUT_MATRIX\OUTPUT_MATRIX);
    fprintf('\ncuff type: %s\n',cuff_types{i})
    fit_and_plot_xy_data(x,y,x_shift,fit_obj,gof,std_type)
    xlabel('frequency (kHz');
    ylabel(thresholds_str);
    title(cuff_types{i})
    
    % change the line color to that specified by line_colors
    gg=get(gca,'Children'); gg(2).Color = line_colors{i}; gg(2).LineWidth=2;
    
end

%%
    %}

    %}
    
end

%% Section 7: Plot raw data

%%% Specify how much to offset the dots so that the electrode type effects
%%% are visible in the raw data
visualization_offset_vals = [-0.25, 0, 0.25];

%%% Specify whether to plot log transformed data or not
flag_plot_log_transformed_data = 0;

%%% For each electrode type, plot the raw data
electrode_type_to_plot = {'monopolar','bipolar','tripolar'};
colors = [
         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    ];


% for onset response data, make a new figure window; otherwise, will just
% superpose data onto the polynomial fits
if (strcmp(metric_str,'onset'))
    custom_figure('position',[ 231         262        1085         418]);
else
    hold on;
end

for i = 1:length(electrode_type_to_plot)
    indices_to_plot = strcmp(frequency_tests_data_table.electrode_type,electrode_type_to_plot{i});
    if (flag_plot_log_transformed_data)
        subset_of_values_to_plot = log_transform_func(metric(indices_to_plot));
    else
        subset_of_values_to_plot = (metric(indices_to_plot));
    end
    plot(visualization_offset_vals(i)+frequency_tests_data_table.frequencies_kHz(indices_to_plot),...
        subset_of_values_to_plot,'.','MarkerSize',20,'color',colors(i,:));
    %     log_transform_func(metric(indices_to_plot)),'.','MarkerSize',20);
    hold on;
end

set(gca,'FontSize',14)
xlabel('frequency (kHz)')
ylabel(thresholds_str)
% title('raw')
% xlim([4 62])
legend(electrode_type_to_plot,'Location','best')


% set xlim to go to zero
xlim_i = get(gca,'XLim');
xlim([0 xlim_i(2)])


% for onset response data, plot with log scale
if (strcmp(metric_str,'onset'))
    % Set scale to be log scale, label ticks, and set ylim bounds
    set(gca,'YScale','log');
    tick_locations = [0.01,0.1,1,10];
    set(gca,'YTick',tick_locations,...
        'YTickLabels',arrayfun(@(x) num2str(x), tick_locations, 'UniformOutput',false))
end

%% Section 8: Plot raw onset responses

%%% Specify ylim manually
warning('manually setting ylim')
ylim_i = [0.0025 10]; % [N*s]

%%% Expand the thresholds table by KHF amplitude relative to block
%%% threshold

if (strcmp(metric_str,'onset'))
    if (~strcmp(average_type,'alt_median'))
        warning('This plotting code is most robust when the alt_median is used so that the onset responses plotted correspond to the interpolated onset responses plotted');
        pause(0.5)
    end
    
    electrode_type_to_plot = {'monopolar','bipolar','tripolar'};
    for electrode_type_ind = 1:length(electrode_type_to_plot)
        
        %%% Specify electrode type to plot
        cuff_type = electrode_type_to_plot{electrode_type_ind};
        
        %%% Specify whether to plot raw or normalized data
        amplitude_type = 'normalized';
        
        %%% Determine variable to use basd on amplitude_type specified
        if strcmp(amplitude_type,'normalized')
            amplitudes_metric = frequency_tests_data_table.normalized_KHF_amplitudes_all;
        elseif strcmp(amplitude_type,'raw')
            amplitudes_metric = frequency_tests_data_table.KHF_amplitudes_all;
        else
            error('invalid value for amplitude_type')
        end
        
        %%% Select unique color map for each nerve of matching cuff type
        nerves_of_matching_cuff_type = unique(frequency_tests_data_table.expdate(find(strcmp(frequency_tests_data_table.electrode_type,cuff_type))));
        assert(length(nerves_of_matching_cuff_type)<=7,'The manually set color map below assumes there are at most seven nerves for any given cuff type');
        color_map_matrix = [
            0    0.4470    0.7410
            0.8500    0.3250    0.0980
            0.9290    0.6940    0.1250
            0.4940    0.1840    0.5560
            0.4660    0.6740    0.1880
            0.3010    0.7450    0.9330
            0.6350    0.0780    0.1840
            ];
        nerve_color_map = containers.Map(nerves_of_matching_cuff_type,num2cell(...
            color_map_matrix(1:length(nerves_of_matching_cuff_type),:),2));
        
        %%% For every frequency, plot all nerves in a panel
        frequencies_all = 10:10:60; % [kHz]
        figure('color',[1 1 1],'position',[685,32,422,964])
        for freq_idx = 1:length(frequencies_all)
            %         subplot(length(frequencies_all),1,freq_idx)
            subaxis(length(frequencies_all),1,freq_idx,'ML',0.2,'MR',0.025,'MT',0.02,'MB',0.13,'SpacingVert',0.004)
            %%% Iterate through all the rows that have the specified electrode type and
            %%% frequency; plot all the nerves on the same data plot
            matching_rows = find(strcmp(frequency_tests_data_table.electrode_type,cuff_type) & frequency_tests_data_table.frequencies_kHz==frequencies_all(freq_idx));
            
            % sort matching_rows by expdate so that the order of the legend is
            % in expdate order
            [~,sort_order] = sort(frequency_tests_data_table.expdate(matching_rows));
            sorted_matching_rows = matching_rows(sort_order);
            
            for i = 1:length(sorted_matching_rows)
                row_idx = sorted_matching_rows(i);
                
                plot(amplitudes_metric{row_idx},frequency_tests_data_table.onset_responses_all{row_idx},'.-','MarkerSize',20,'LineWidth',2,...
                    'Color',nerve_color_map(frequency_tests_data_table.expdate{row_idx}));
                hold on;
            end
            
            %         xlim_i = get(gca,'XLim');
            %         xlim([0.9 max(xlim_i)])
            
            % Set a common x axis bounds, and only label the bottom panel's x axis
            xlim([0.9 3.5])
            if (freq_idx==length(frequencies_all))
                xlabel('KHF amplitude (divided by block thresh)')
            else
                set(gca,'XTickLabels','')
            end
            
            ylabel({sprintf('%d kHz',frequencies_all(freq_idx)),'onset resp (N*s)'})
            if (freq_idx==1)
                title(cuff_type)
            end
            set(gca,'FontSize',11)
            
            % Set scale to be log scale, label ticks, and set ylim bounds
            set(gca,'YScale','log');
            tick_locations = [0.01,0.1,1,10];
            set(gca,'YTick',tick_locations,...
                'YTickLabels',arrayfun(@(x) num2str(x), tick_locations, 'UniformOutput',false))
            ylim(ylim_i)
        end
        
        warning('code assumes the plots are in expdate order; legend will label accordingly')
        legend_labels = unique(frequency_tests_data_table.expdate(strcmp(frequency_tests_data_table.electrode_type,cuff_type)));
        legend(legend_labels,'NumColumns',3);
    end
end


%% Section 9: For each nerve, get & plot block threshold at minimum frequency vs. block threshold at frequency minimum block threshold

%%% Identify unique nerves
unique_nerves = unique(frequency_tests_data_table.expdate);

% Ignore the nerves that are marked for exclusion from the min freq
% analysis
metadata_table = readtable('minimum_blocking_frequency_posthoc_analysis_v3.xlsx');
nerve_to_inclusion_map = containers.Map(metadata_table.ExperimentDate,metadata_table.IncludeMinFreqAnalysis);
unique_nerves = unique_nerves(find(cellfun(@(x) nerve_to_inclusion_map(x), unique_nerves)));

%%% For each nerve, find the minimum frequency and its energy and onset
%%% response
min_freq_data_table = table(zeros(length(unique_nerves),1),zeros(length(unique_nerves),1),zeros(length(unique_nerves),1),...
    'VariableNames',{'frequency_kHz','power_mW','onset_response_Ns'});
for i = 1:length(unique_nerves)
    row_indices = find(strcmp(frequency_tests_data_table.expdate,unique_nerves{i}));
    
    % Offset the frequencies differently for each cuff type for subsequent
    % visualization
    %%% The first element should has the same electrode type as all others,
    %%% so extract the electrode type based on the first element
    electrode_type_i = frequency_tests_data_table.electrode_type{row_indices(1)};
    if (strcmp(electrode_type_i,'monopolar'))
        offset_freq_for_visualization_Hz = +0.05; % [Hz]
    elseif (strcmp(electrode_type_i,'bipolar'))
        offset_freq_for_visualization_Hz = 0; % [Hz]
    elseif (strcmp(electrode_type_i,'tripolar'))
        offset_freq_for_visualization_Hz = -0.05; % [Hz]
    else
        error('unexpected value of electrode_type_i: %s',electrode_type_i);
    end
    [min_val,idx] = min(frequency_tests_data_table.frequencies_kHz(row_indices));
    min_freq_data_table.frequency_kHz(i) = offset_freq_for_visualization_Hz + min_val;
    
    min_freq_data_table.power_mW(i) = frequency_tests_data_table.thresholds_amplitude_mW(row_indices(idx));
    min_freq_data_table.onset_response_Ns(i) = frequency_tests_data_table.onset_response_at_single_amp_Ns(row_indices(idx));
end

%%% For each nerve, find the minimum energy and the frequency and onset
%%% response at the minimum energy
min_energy_data_table = table(zeros(length(unique_nerves),1),zeros(length(unique_nerves),1),zeros(length(unique_nerves),1),...
    'VariableNames',{'frequency_kHz','power_mW','onset_response_Ns'});
for i = 1:length(unique_nerves)
    row_indices = find(strcmp(frequency_tests_data_table.expdate,unique_nerves{i}));
    [min_energy_data_table.power_mW(i),idx] = min(frequency_tests_data_table.thresholds_amplitude_mW(row_indices));
    
    % Offset the frequencies differently for each cuff type for subsequent
    % visualization
    %%% The first element should has the same electrode type as all others,
    %%% so extract the electrode type based on the first element
    electrode_type_i = frequency_tests_data_table.electrode_type{row_indices(1)};
    if (strcmp(electrode_type_i,'monopolar'))
        offset_freq_for_visualization_Hz = +0.05; % [Hz]
    elseif (strcmp(electrode_type_i,'bipolar'))
        offset_freq_for_visualization_Hz = 0; % [Hz]
    elseif (strcmp(electrode_type_i,'tripolar'))
        offset_freq_for_visualization_Hz = -0.05; % [Hz]
    else
        error('unexpected value of electrode_type_i: %s',electrode_type_i);
    end
    min_energy_data_table.frequency_kHz(i) = offset_freq_for_visualization_Hz + frequency_tests_data_table.frequencies_kHz(row_indices(idx));
    
    min_energy_data_table.onset_response_Ns(i) = frequency_tests_data_table.onset_response_at_single_amp_Ns(row_indices(idx));
end

%%% For each nerve, find the minimum onset response and the frequency and
%%% energy at the minimum onset response
min_onset_data_table = table(zeros(length(unique_nerves),1),zeros(length(unique_nerves),1),zeros(length(unique_nerves),1),...
    'VariableNames',{'frequency_kHz','power_mW','onset_response_Ns'});
for i = 1:length(unique_nerves)
    row_indices = find(strcmp(frequency_tests_data_table.expdate,unique_nerves{i}));
    [min_onset_data_table.onset_response_Ns(i),idx] = min(frequency_tests_data_table.onset_response_at_single_amp_Ns(row_indices));
    
    % Offset the frequencies differently for each cuff type for subsequent
    % visualization
    %%% The first element should has the same electrode type as all others,
    %%% so extract the electrode type based on the first element
    electrode_type_i = frequency_tests_data_table.electrode_type{row_indices(1)};
    if (strcmp(electrode_type_i,'monopolar'))
        offset_freq_for_visualization_Hz = +0.05; % [Hz]
    elseif (strcmp(electrode_type_i,'bipolar'))
        offset_freq_for_visualization_Hz = 0; % [Hz]
    elseif (strcmp(electrode_type_i,'tripolar'))
        offset_freq_for_visualization_Hz = -0.05; % [Hz]
    else
        error('unexpected value of electrode_type_i: %s',electrode_type_i);
    end
    min_onset_data_table.frequency_kHz(i) = offset_freq_for_visualization_Hz + frequency_tests_data_table.frequencies_kHz(row_indices(idx));
    
    min_onset_data_table.power_mW(i) = frequency_tests_data_table.thresholds_amplitude_mW(row_indices(idx));
end

%%% Get electrode identity to group the plots by electrode type
electrode_types = cellfun(@(x) ...
    frequency_tests_data_table.electrode_type{find(strcmp(frequency_tests_data_table.expdate,x))}, unique_nerves,...
    'UniformOutput',false);

%%% Plot the various combos
xy_combos = {
    min_freq_data_table.frequency_kHz, min_energy_data_table.frequency_kHz,'minimum blocking freq (kHz)','most efficient blocking frequency (kHz)'
    min_freq_data_table.power_mW, min_energy_data_table.power_mW,'thresh50 (mW) at minimum blocking freq','thresh50 (mW) at most efficient blocking freq'
    min_freq_data_table.onset_response_Ns, min_energy_data_table.onset_response_Ns,'onset resp (N*s) at min blocking freq','onset resp (N*s) at most efficient freq'
    min_freq_data_table.frequency_kHz, min_onset_data_table.frequency_kHz,'freq (kHz) at min blocking freq','freq (kHz) at lowest onset freq'
    min_freq_data_table.power_mW, min_onset_data_table.power_mW,'thresh (mW) at min blocking freq','thresh (mW) at lowest onset freq'
    min_freq_data_table.onset_response_Ns, min_onset_data_table.onset_response_Ns,'onset resp (N*s) at min blocking freq','onset resp (N*s) at lowest onset freq'
    min_energy_data_table.frequency_kHz, min_onset_data_table.frequency_kHz,'freq (kHz) at most efficient freq','freq (kHz) at lowest onset freq'
    min_energy_data_table.power_mW, min_onset_data_table.power_mW,'thresh (mW) at most efficient freq','thresh (mW) at lowest onset freq'
    min_energy_data_table.onset_response_Ns, min_onset_data_table.onset_response_Ns,'onset resp (N*s) at most efficient freq','onset resp (N*s) at lowest onset freq'
};

for i = 1:size(xy_combos,1)
    x = xy_combos{i,1};
    y = xy_combos{i,2};
    xlabel_str = xy_combos{i,3};
    ylabel_str = xy_combos{i,4};
    
    custom_figure;
    group_labels = {'monopolar','bipolar','tripolar'};
    group_colors = {[     0    0.4470    0.7410],[0.8500    0.3250    0.0980],[0.9290    0.6940    0.1250]};
    group_data = cellfun(@(group_i_label) ...
        [x(strcmp(electrode_types,group_i_label)),y(strcmp(electrode_types,group_i_label))],...
        group_labels,'UniformOutput',false);
    scatter_pairwise_comparison(group_data,group_colors,xlabel_str,ylabel_str,0,group_labels)
    legend off 
    
    %%% Perform stat tests
    p = signrank(x,y);
    fprintf('%s vs. %s: p = %0.6e\n\n',ylabel_str,xlabel_str,p);
   

end

%% Write the relevant text with auto-filled stats

range_and_median_str_func = @(x) sprintf('(%d to %d kHz; median: %d kHz)', round(min(x)),round(max(x)),round(median(x)));
paired_comparison_a_over_b_func = @(a,b) sprintf('%0.3f [%0.3f, %0.3f]; p=%0.1e', median(a./b),min(a./b),max(a./b),signrank(a,b));
% 0.78 [0.6, 1]; p=4.9e-5
% (5 to 9 kHz; median: 7 kHz)
fprintf(['We compared the frequency, threshold, and onset response for the most efficient blocking frequency and for the minimum blocking frequency. The minimum blocking frequency ', range_and_median_str_func(min_freq_data_table.frequency_kHz),' was usually lower than the frequency at which minimum block threshold occurred ',range_and_median_str_func(min_energy_data_table.frequency_kHz),' ([***Figure 3A***]; median [min, max] of ratios & signed rank test: ',paired_comparison_a_over_b_func(min_freq_data_table.frequency_kHz,min_energy_data_table.frequency_kHz),'). The block thresholds at the most efficient frequencies were smaller than the block thresholds of the minimum blocking frequency ([***Figure 3B***]; ',paired_comparison_a_over_b_func(min_energy_data_table.power_mW,min_freq_data_table.power_mW),'). Onset responses at the most efficient frequencies were comparable to those of the minimum blocking frequencies ([***Figure 3C***]; ',paired_comparison_a_over_b_func(min_energy_data_table.onset_response_Ns,min_freq_data_table.onset_response_Ns),'), although the minimum blocking frequency often had larger onset responses, consistent with prior literature showing higher onset response at lower frequencies [***CITE***]. \n']);
fprintf(['We also compared the frequency, threshold, and onset response at the frequency that generated the smallest onset response to the minimum blocking frequency and to the most efficient blocking frequency. Minimum onset response occurred at a higher frequency ',range_and_median_str_func(min_onset_data_table.frequency_kHz),' than both the minimum blocking frequency (***Figure 3D***; ',paired_comparison_a_over_b_func(min_onset_data_table.frequency_kHz,min_freq_data_table.frequency_kHz),') and the most efficient blocking frequency (***Figure 3G***; ',paired_comparison_a_over_b_func(min_onset_data_table.frequency_kHz,min_energy_data_table.frequency_kHz),'). At those higher frequencies, onset response was substantially smaller compared to the minimum blocking frequency (***Figure 3F***; ',paired_comparison_a_over_b_func(min_onset_data_table.onset_response_Ns,min_freq_data_table.onset_response_Ns),') and compared to the most efficient blocking frequency (***Figure 3I***; ',paired_comparison_a_over_b_func(min_onset_data_table.onset_response_Ns,min_energy_data_table.onset_response_Ns),'). Meanwhile, block thresholds that generated the minimum onset response were much larger compared to the minimum blocking frequency (***Figure 3E***; ',paired_comparison_a_over_b_func(min_onset_data_table.power_mW,min_freq_data_table.power_mW),') and compared to the most efficient blocking frequency (***Figure 3H***; ',paired_comparison_a_over_b_func(min_onset_data_table.power_mW,min_energy_data_table.power_mW),'). Results were similar when using Method 2 to define minimum blocking frequency (***Supplementary Figure 2***). \n']);

%% Linear models processing
% If treating frequency as categorical, only include a data point if
% there are at least three data points of the same frequency and at least
% two different cuff types
if (strcmp(regression_parameters.frequency_type,'categorical'))
    unique_frequencies = unique(frequency_tests_data_table.frequencies_kHz);
    for freq_idx = 1:length(unique_frequencies)
        % check whether there are at least three data points at this frequency
        n_points = sum(frequency_tests_data_table.frequencies_kHz==unique_frequencies(freq_idx));
        is_n_points_too_small = n_points<3;
        
        % check whether there are at least two different cuff types at this
        % frequency
        n_cuff_types = length(unique(frequency_tests_data_table.electrode_type(frequency_tests_data_table.frequencies_kHz==unique_frequencies(freq_idx))));
        is_n_cuff_types_too_small = n_cuff_types<2;
        
        % check whether there are at least two different nerves at this
        % frequency
        n_nerves = length(unique(frequency_tests_data_table.expdate(frequency_tests_data_table.frequencies_kHz==unique_frequencies(freq_idx))));
        is_n_nerves_too_small = n_nerves<2;
        
        assert(n_nerves==n_points,'Since removing replicates, n_nerves shoudl always equal n_points')
        
        % if either of the above conditions failed, exclude all data points at
        % this frequency
        if (is_n_points_too_small || is_n_cuff_types_too_small || is_n_nerves_too_small)
            warning('removing %d kHz because not enough data points or not enough cuff types or not enough nerves there',unique_frequencies(freq_idx))
            rows_to_remove = find(frequency_tests_data_table.frequencies_kHz==unique_frequencies(freq_idx));
            frequency_tests_data_table(rows_to_remove,:) = [];
            metric(rows_to_remove) = [];
        end
    end
end


% % Perform regression

% set up regression input
expdate_indices = cellfun(@(x) find(strcmp(x,unique(frequency_tests_data_table.expdate))), frequency_tests_data_table.expdate);
expdate_variable_matrix = [zeros(1,max(expdate_indices)-1); eye(max(expdate_indices)-1)];
regression_input_nerve = expdate_variable_matrix(expdate_indices,:);
regression_input_nerve_str = arrayfun(@(x) num2str(x,'exp%d'),2:(length(unique(frequency_tests_data_table.expdate))),'UniformOutput',false);

electrode_type_indices = cellfun(@(x) find(strcmp(x,unique(frequency_tests_data_table.electrode_type))), frequency_tests_data_table.electrode_type);
electrode_type_variable_matrix = [zeros(1,max(electrode_type_indices)-1); eye(max(electrode_type_indices)-1)];
regression_input_elect = electrode_type_variable_matrix(electrode_type_indices,:);
% warning('attempting to get e1 too.. (even though %TODO the two lines below are the only edits I''ve made toward this')
% electrode_type_variable_matrix = [eye(max(electrode_type_indices))];
% regression_input_elect = electrode_type_variable_matrix(electrode_type_indices,:);
regression_input_elect_str = arrayfun(@(x) num2str(x,'elec%d'),2:(length(unique(frequency_tests_data_table.electrode_type))),'UniformOutput',false);

if (strcmp(regression_parameters.frequency_type,'categorical'))
    
    % Specify basis frequency
    basis_freq = 10; % [kHz]
    
    % Compute unique frequencies
    unique_frequencies = unique(frequency_tests_data_table.frequencies_kHz);
    
    % Resort the unique frequency vector so that the desired basis
    % frequency is the first element in the vector; this will ensure that
    % the lines that follow assign a zero vector (i.e., no coefficients) to
    % the linear model when that basis frequency is chosen
    idx = find(basis_freq==unique_frequencies);
    assert(~isempty(idx),'Specified basis_freq must be one of the frequencies sampled');
    resorted_unique_frequencies = [unique_frequencies(idx);setdiff(unique_frequencies,unique_frequencies(idx))];
    
    frequency_indices = arrayfun(@(x) find(x==resorted_unique_frequencies), frequency_tests_data_table.frequencies_kHz);
    frequency_variable_matrix = [zeros(1,max(frequency_indices)-1); eye(max(frequency_indices)-1)];
    regression_input_freq = frequency_variable_matrix(frequency_indices,:);
    regression_input_freq_str = arrayfun(@(x) num2str(x,'freq%d'),2:(length(resorted_unique_frequencies)),'UniformOutput',false);
    
    % elseif (strcmp(regression_parameters.frequency_type,'categorical_with_interactions'))
    
else
    split_frequencies = 9; % [kHz]
    
    is_left_freq = frequency_tests_data_table.frequencies_kHz<=split_frequencies(1);
    is_right_freq = frequency_tests_data_table.frequencies_kHz>split_frequencies(1);
    
    if (regression_parameters.flag_log_transform)
        tmp_freq = log_transform_func(frequency_tests_data_table.frequencies_kHz);
    else
        tmp_freq = frequency_tests_data_table.frequencies_kHz;
    end
    regression_input_freq = ([tmp_freq.*is_left_freq, tmp_freq.*is_right_freq]);
    regression_input_freq_str = {'freq_left','freq_right'};
    
    warning('for continous freq, manually overriding split freq and using instead a single freq range')
    pause(1)
    regression_input_freq = tmp_freq; %frequency_tests_data_table.frequencies_kHz;
    regression_input_freq_str = {'freq'};
    
    %{
    warning('for continous freq, manually overriding split freq and using instead a quadratic function of freq range')
    pause(1)
    regression_input_freq = [frequency_tests_data_table.frequencies_kHz,frequency_tests_data_table.frequencies_kHz.^2];
    regression_input_freq_str = {'freq','freq_squared'};
    %}
end

if (regression_parameters.include_electrode_effects && ~regression_parameters.include_random_effects)
    regression_input = [regression_input_freq,regression_input_elect];
    regression_input_str = [...
        regression_input_freq_str,...
        regression_input_elect_str];
elseif (regression_parameters.include_random_effects && ~regression_parameters.include_electrode_effects)
    regression_input = [regression_input_nerve,regression_input_freq];
    regression_input_str = [regression_input_nerve_str,...
        regression_input_freq_str];
    columns_to_keep = 1:length(regression_input_nerve_str);
else
    
    % since doing both nerve and electrode effects, need to remove one of
    % the nerve effects for each electrode type, since cannot adjust for
    % *every* nerve if adjusting for electrode type also
    monopolar_idx = 1; % monopolar
    tripolar_idx = 3; % tripolar
    columns_to_remove = [...
        find(expdate_variable_matrix(expdate_indices(find(electrode_type_indices==monopolar_idx,1,'first')),:)),...
        find(expdate_variable_matrix(expdate_indices(find(electrode_type_indices==tripolar_idx,1,'first')),:))];
    columns_to_keep = setdiff(1:length(regression_input_nerve_str),columns_to_remove);
    
    regression_input = [regression_input_nerve(:,columns_to_keep),regression_input_freq,regression_input_elect];
    regression_input_str = [regression_input_nerve_str(columns_to_keep),...
        regression_input_freq_str,regression_input_elect_str];
end

% include interaction effects
flag_include_continuous_interactions = 0;
if (flag_include_continuous_interactions)
    warning('including interaction effects')
    assert(strcmp(regression_parameters.frequency_type,'continuous'),'interaction effects only work with continuous freq right now')
    regression_input_interaction = (electrode_type_indices==3).*frequency_tests_data_table.frequencies_kHz;
    regression_input = [regression_input, regression_input_interaction];
    regression_input_interaction_str = {'freq_interaction'};
    regression_input_str = [regression_input_str,regression_input_interaction_str];
end

if (flag_include_hybrid_interactions)
    assert(strcmp(regression_parameters.frequency_type,'categorical'),'must use categorical freq')

    warning('manually set cutoff freq to %d kHz by manually sweeping through freq values from 9 to 29 in 1 kHz increments and picking the freq that maximized Adjusted R-Squared',freq_cutoff);
    regression_input_interaction = (electrode_type_indices==3) & (frequency_tests_data_table.frequencies_kHz>freq_cutoff);
    regression_input = [regression_input, regression_input_interaction];
    regression_input_interaction_str = {'freq_interaction'};
    regression_input_str = [regression_input_str,regression_input_interaction_str];
    % add an interaction term that is only present at a certain subset of
    % frequencies above a certain kilohertz value
    % in v1, specify the freqeuency cutoff; in v2, find the optimal cutoff
    % empirically
    
    %%%
end

% set up regression output
regression_output = metric;
regression_output_str = {'threshold'};

if (regression_parameters.flag_log_transform)
    regression_output = log_transform_func(regression_output);
    regression_output_str{1} = ['logX_',regression_output_str{1}];
    regression_input_str = cellfun(@(x) ['logX_',x], regression_input_str, 'UniformOutput', false);
end

%{
        same as:
            [ones(size(regression_input,1),1),regression_input]\regression_output
        but automatically computes stats: [R^2,F,p,SSE]
%}

if (regression_parameters.weighting)
    local_reweight_loop_control_flag = true;
    coefficients_manual = [];
    reweight_count = 0;
    while (local_reweight_loop_control_flag)
        % if reweighting is not specified, then just run the below once
        if (~regression_parameters.reweight)
            local_reweight_loop_control_flag = false;
        end
        
        %%% If there reweighting is specified and one iteration of the loop
        %%% already completed...
        if (~isempty(coefficients_manual) && regression_parameters.reweight)
            uncertainty = ([ones(size(regression_input,1),1),regression_input]*coefficients_manual).^1;
            %             sigma = 1./(inv_log_transform_func([ones(size(regression_input,1),1),regression_input]*coefficients_manual));
            reweight_count = reweight_count + 1;
            if (reweight_count==reweight_count_max)
                local_reweight_loop_control_flag = false;
            end
            
            %%% Otherwise, set initial weights for first iteration (i.e.,
            %%% the only iteration if reweighting flag is false)
        else
            flag_init_reweighting_scheme = 2;
            switch flag_init_reweighting_scheme
                %%% Use the residual of a simple linear fit (freq vs.
                %%% output) to define uncertainty values
                case 1
                    % make sigma (i.e., expected standard deviation of the data
                    % point) be greater as the value of the data goes up by setting
                    % based on a simple fit
                    x = regression_input_freq(:,1);
                    y = regression_output;
                    simple_linear_fit = [x, ones(size(x))]\y;
                    uncertainty = std((polyval(simple_linear_fit,x)).^1 - y);
                    %{
                        figure; plot(x,sigma)
                        hold on; plot(regression_input_freq(:,1),regression_output,'.')
                    %}
                    
                    %             sigma = 1./(inv_log_transform_func(regression_output)); % make sigma (i.e., expected standard deviation of the data point) be greater as the value of the data goes up
                    
                case 2
                    %%% Calculate the variance of the data at each frequency
                    %%% and electrode type combination for frequencies 10, 20,
                    %%% 30, 40, 50, and 60 kHz, then we fit a line to the
                    %%% variance for each electrode type and interpolated or
                    %%% extrapolated variance values at other frequencies
                    %%% Do this separately for each electrode type
                    uncertainty = ones(size(regression_output));
                    
%                     warning('disabled weighting...'); pause(2)
                    
                    for cuff_type_idx = 1:length(cuff_types)
                        rows_to_fill = find(strcmp(cuff_types{cuff_type_idx},frequency_tests_data_table.electrode_type));
                        training_freq = 10:10:60; % [kHz]
                        training_std = zeros(size(training_freq));
                        for freq_idx = 1:length(training_freq)
                            training_data_rows = find(ismember(frequency_tests_data_table.frequencies_kHz,training_freq(freq_idx)) & ...
                                strcmp(cuff_types{cuff_type_idx},frequency_tests_data_table.electrode_type));
                            training_std(freq_idx) = std(regression_output(training_data_rows));
                        end
                        
                        %%% Interpolate value (using 10 kHz for frequencies
                        %%% less than 10 kHz so that the code below
                        %%% effectively performed nearest neighbors
                        %%% interpolation and nearest neighbors
                        %%% extrapolation
                        uncertainty(rows_to_fill) = interp1(training_freq,training_std,max(frequency_tests_data_table.frequencies_kHz(rows_to_fill),10),...
                            'nearest',NaN);
                        assert(~any(isnan(uncertainty)),'sigma should not have any NaN values at this point')
                        
                        %%% Define the weight as the inverse of the
                        %%% uncertainty (since a more certain value 
                    end
                    
                    
                otherwise
                    error('Invalid value for flag_init_reweighting_scheme')
            end
        end
        
        % modify regression_input and regression_output using the weights
        regression_input_final = [(ones(size(regression_input,1),1)./uncertainty),(regression_input./uncertainty)]; % append a 'ones' column since using fitlm but still doing weighting, so wil set to 'Intercept' false but still having a manually inseted intercept column
        regression_output_final = regression_output./uncertainty;
        
        % set up regression table with variable names
        regression_table = array2table([regression_input_final,regression_output_final],'VariableNames',...
            ['intercept',regression_input_str,regression_output_str]);
        
        coefficients_manual=regression_input_final\regression_output_final;
        
        % perform regression
%         SSE = sum(([ones(size(regression_input,1),1),regression_input]*coefficients_manual - regression_output).^2);
        mdl = fitlm(regression_table,'Intercept',false); % same as (regression_input_final\regression_output_final)'

    end
    
    %%% Calculate residuals
    flag_use_weighted_data_for_residuals = 1;
    if (flag_use_weighted_data_for_residuals)
        residuals = regression_input_final*mdl.Coefficients.Estimate - regression_output_final;
    else
        residuals = [ones(size(regression_input,1),1),regression_input]*mdl.Coefficients.Estimate - regression_output;
    end
else
    regression_input_final = regression_input;
    regression_output_final = regression_output;

    % set up regression table with variable names
    regression_table = array2table([regression_input_final,regression_output_final],'VariableNames',...
        [regression_input_str,regression_output_str]);
    
    % perform regression
    mdl = fitlm(regression_table);
    
    %%% Calculate residuals
    residuals = [ones(size(regression_input_final,1),1),regression_input_final]*mdl.Coefficients.Estimate - regression_output;
end

%%% Display the linear model
disp(mdl)

%%% Plot residual histogram and run Anderson-Darling test as gut check for the normality of the fit
figure('color',[1 1 1]); histogram(residuals); xlabel('residual'); ylabel('count'); set(gca,'FontSize',14)
fprintf('Anderson-Darling test output (null hypothesis is that residuals are normally distributed):')
[h,p,~,~]=adtest(residuals)

%%% Extract the coefficients and confidence intervals
confidence_intervals = mdl.coefCI;

n = 1; % index at which to start parsing for mdl outputs
intercept = mdl.Coefficients.Estimate(n);
intercept_CI = confidence_intervals(n,:);

if (regression_parameters.include_random_effects)
    n = n(end)+(1:length(regression_input_nerve_str(columns_to_keep)));
    nerve_coeffs = mdl.Coefficients.Estimate(n);
    nerve_coeffs_CI = confidence_intervals(n,:);
    adjustment_for_nerve = regression_input_nerve(:,columns_to_keep)*nerve_coeffs;
else
    adjustment_for_nerve = zeros(size(regression_input,1),1);
end

n = n(end)+(1:length(regression_input_freq_str));
frequency_coeffs = mdl.Coefficients.Estimate(n);
frequency_coeffs_CI = confidence_intervals(n,:);
adjustment_for_freq = regression_input_freq*frequency_coeffs;

freq_readjustment_factor = 0;


if (regression_parameters.include_electrode_effects)
    n = n(end)+(1:length(regression_input_elect_str));
    elect_coeffs = mdl.Coefficients.Estimate(n);
    elect_coeffs_CI = confidence_intervals(n,:);
    
    adjustment_for_elect = regression_input_elect*elect_coeffs;
else
    adjustment_for_elect = zeros(size(regression_input,1),1);
end


if (flag_include_continuous_interactions)
    n = n(end)+(1:length(regression_input_interaction_str));
    interaction_term = mdl.Coefficients.Estimate(n);
    interaction_term_CI = confidence_intervals(n,:);
    adjustment_for_interaction_term = regression_input_interaction*interaction_term;
    %     interaction_term_readjustment = freq_for_readjustment*interaction_term.*(regression_input_interaction~=0); % by default, adjusting for the interaction term adjusts the threshold to what they would theoretically be at 0 kHz, so re-adjust for 10 kHz instead
elseif (flag_include_hybrid_interactions)
    n = n(end)+(1:length(regression_input_interaction_str));
    interaction_term = mdl.Coefficients.Estimate(n);
    interaction_term_CI = confidence_intervals(n,:);
    adjustment_for_interaction_term = regression_input_interaction*interaction_term;
else
    adjustment_for_interaction_term = 0;
end

output_adjusted_for_nerve_and_elect = regression_output-adjustment_for_nerve-adjustment_for_elect-adjustment_for_interaction_term...
    +freq_readjustment_factor;
if (regression_parameters.flag_log_transform)
    output_adjusted_for_nerve_and_elect = inv_log_transform_func(output_adjusted_for_nerve_and_elect);
end


output_adjusted_for_nerve_only = regression_output-adjustment_for_nerve...
    +freq_readjustment_factor;
if (regression_parameters.flag_log_transform)
    output_adjusted_for_nerve_only = inv_log_transform_func(output_adjusted_for_nerve_only);
end


% Adjustment for nerve and freq, but *not* for interaction terms because
% the subsequent code displays the interaction effects graphically, so if
% they get removed here they will not show up
output_adjusted_for_nerve_and_freq = regression_output-adjustment_for_nerve-adjustment_for_freq...
    +freq_readjustment_factor; %+interaction_term_readjustment
if (regression_parameters.flag_log_transform)
    output_adjusted_for_nerve_and_freq = inv_log_transform_func(output_adjusted_for_nerve_and_freq);
end

% Adjustment for nerve and freq and interaction terms 
output_adjusted_for_nerve_and_freq_and_interaction = regression_output-adjustment_for_nerve-adjustment_for_freq...
    +freq_readjustment_factor-adjustment_for_interaction_term; %+interaction_term_readjustment
if (regression_parameters.flag_log_transform)
    output_adjusted_for_nerve_and_freq_and_interaction = inv_log_transform_func(output_adjusted_for_nerve_and_freq_and_interaction);
end


%%% If applicable, print out electrode coefficients and confidence intervals
if (regression_parameters.include_electrode_effects)
    fprintf('electrode coefficients and confidence intervals:\n')
    if (regression_parameters.flag_log_transform)
        disp(inv_log_transform_func(elect_coeffs))
        disp(inv_log_transform_func(elect_coeffs_CI))
        
        if (flag_include_continuous_interactions || flag_include_hybrid_interactions)
            disp(inv_log_transform_func(interaction_term))
            disp(inv_log_transform_func(interaction_term_CI))
        end
    else
        disp(elect_coeffs)
        disp(elect_coeffs_CI)
        
        if (flag_include_continuous_interactions || flag_include_hybrid_interactions)
            disp(interaction_term)
            disp(interaction_term_CI)
        end
    end
end


%% Plot data after adjusting for nerve and electrode effects

%%% Specify how much to offset the dots so that the electrode type effects
%%% are visible in the raw data
visualization_offset_vals = [-0.25, 0, 0.25];

%%% Specify whether to plot log transformed data or not
flag_plot_log_transformed_data = 0;

%%% Specify which dataset to plot
flag_dataset_to_plot = 2;
switch flag_dataset_to_plot
    case 0
        values_to_plot = metric;
        value_str = 'raw (subset used for linear model)';
    case 1
        values_to_plot = output_adjusted_for_nerve_only;
        value_str = 'adjusted for nerve only';
    case 2
        values_to_plot = output_adjusted_for_nerve_and_elect;
        value_str = 'adjusted for nerve and electrode';
    case 3
        values_to_plot = output_adjusted_for_nerve_and_freq;
        value_str = 'adjusted for nerve and frequency';
    case 4
        values_to_plot = output_adjusted_for_nerve_and_freq_and_interaction;
        value_str = 'adjusted for nerve and frequency';
    otherwise
        error('Invalid value for flag_dataset_to_plot')
end

%%% For each electrode type, plot the raw data
electrode_type_to_plot = {'monopolar','bipolar','tripolar'};
custom_figure('position',[ 231         552        1085         418]);
for i = 1:length(electrode_type_to_plot)
    indices_to_plot = strcmp(frequency_tests_data_table.electrode_type,electrode_type_to_plot{i});
    if (flag_plot_log_transformed_data)
        subset_of_values_to_plot = log_transform_func(values_to_plot(indices_to_plot));
        ylabel_str_prefix = sprintf('log_%d of ',regression_parameters.log_base);
    else
        subset_of_values_to_plot = (values_to_plot(indices_to_plot));
        ylabel_str_prefix = '';
    end
    plot(visualization_offset_vals(i)+frequency_tests_data_table.frequencies_kHz(indices_to_plot),...
        subset_of_values_to_plot,'.','MarkerSize',20);
    %     log_transform_func(metric(indices_to_plot)),'.','MarkerSize',20);
    hold on;
end

set(gca,'FontSize',14)
xlabel('frequency (kHz)')
ylabel([ylabel_str_prefix, ' ' thresholds_str])
title(value_str)
xlim([4 62])
legend(electrode_type_to_plot,'Location','North')


%% Plot data after adjusting for nerve and frequency effects
fig_handle = custom_figure('position',[230,41,917,418]);

colors_cell_array = {
    frequency_tests_data_table.frequencies_kHz(strcmp(frequency_tests_data_table.electrode_type,'monopolar'))
    frequency_tests_data_table.frequencies_kHz(strcmp(frequency_tests_data_table.electrode_type,'bipolar'))
    frequency_tests_data_table.frequencies_kHz(strcmp(frequency_tests_data_table.electrode_type,'tripolar') & ...
    frequency_tests_data_table.frequencies_kHz<=freq_cutoff)
    frequency_tests_data_table.frequencies_kHz(strcmp(frequency_tests_data_table.electrode_type,'tripolar') & ...
    frequency_tests_data_table.frequencies_kHz>freq_cutoff)
    };

% Add compare_group_fitness to path
scatter_size = 30;
compare_group_fitness({...
    output_adjusted_for_nerve_and_freq(strcmp(frequency_tests_data_table.electrode_type,'monopolar')),...
    output_adjusted_for_nerve_and_freq(strcmp(frequency_tests_data_table.electrode_type,'bipolar')),...
    output_adjusted_for_nerve_and_freq(strcmp(frequency_tests_data_table.electrode_type,'tripolar') & ...
    frequency_tests_data_table.frequencies_kHz<=freq_cutoff),...
    output_adjusted_for_nerve_and_freq(strcmp(frequency_tests_data_table.electrode_type,'tripolar') & ...
    frequency_tests_data_table.frequencies_kHz>freq_cutoff)...
    },fig_handle,1,colors_cell_array,scatter_size);
h = colorbar;
ylabel(h,'frequency (kHz)')

set(gca,'FontSize',14)
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'monopolar','bipolar',['tripolar <=', sprintf('%dk',freq_cutoff)],['tripolar >',sprintf('%dk',freq_cutoff)]},'XTickLabelRotation',0)
ylabel(thresholds_str)
xlim([0 5])
fprintf('regression complete\n');


