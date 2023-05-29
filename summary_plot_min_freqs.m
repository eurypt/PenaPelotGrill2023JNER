%{

Plot Figure 2 of manuscript comparing minimum frequency across different
methods and across different electrode types.

%}

function summary_plot_min_freqs()

% add path to internal dependencies (i.e., functions I wrote to facilitate
% plotting/analysis/etc.)
addpath('internal_dependencies');

rng('default')

% Load manually determined minimum blocking frequencies
metadata_table = readtable('minimum_blocking_frequency_posthoc_analysis_v3.xlsx');

% Figure 2A
% Specify whether or not to use manually determined minimum frequencies
flag_use_manual_min_thresholds = 1;
% Specify which automatic method to use for determining minimum frequency;
% set to either 'max' or 'AUC'; only use if
% flag_use_manual_min_thresholds is set to 0
auto_method = [];
% Make comparison plot
make_summary_plot(metadata_table,flag_use_manual_min_thresholds,auto_method)

% Figure 2B
% Make additional comparisons
flag_use_manual_min_thresholds = 0;
auto_method = 'AUC';
make_summary_plot(metadata_table,flag_use_manual_min_thresholds,auto_method);

% Figure 2C
flag_use_manual_min_thresholds = 0;
auto_method = 'max';
make_summary_plot(metadata_table,flag_use_manual_min_thresholds,auto_method);

function make_summary_plot(metadata_table,flag_use_manual_min_thresholds,auto_method)
if (flag_use_manual_min_thresholds)
    % If the data point is marked for inclusion in the min freq
    % analysis, then extract the min freq
    min_block_thresholds_monopolar = metadata_table.MinimumBlockingFreqKilohertz(strcmp(metadata_table.CuffType,'Monopolar') & metadata_table.IncludeMinFreqAnalysis);
    min_block_thresholds_bipolar = metadata_table.MinimumBlockingFreqKilohertz(strcmp(metadata_table.CuffType,'Bipolar') & metadata_table.IncludeMinFreqAnalysis);
    min_block_thresholds_tripolar = metadata_table.MinimumBlockingFreqKilohertz(strcmp(metadata_table.CuffType,'Tripolar') & metadata_table.IncludeMinFreqAnalysis);
else
    load(['frequency_tests_data_table_',auto_method,'.mat'],'frequency_tests_data_table')
    for i = 1:size(frequency_tests_data_table,1)
        if (frequency_tests_data_table.frequencies_kHz(i)==1)
            frequency_tests_data_table.thresholds_amplitude_mA(i) = NaN;  % set thresholds to NaN whereever 1 kHz was tested, since none of these blocked, but some spuriously show up as have a valid threshold
        end
    end
    min_block_thresholds_monopolar = [];
    min_block_thresholds_bipolar = [];
    min_block_thresholds_tripolar = [];
    for i = 1:size(metadata_table,1)
        % If the data point is marked for inclusion in the min freq
        % analysis, then extract the min freq 
        if (metadata_table.IncludeMinFreqAnalysis(i)) 
            threshold_rows = find(cellfun(@(x) strcmp(x(1:8),metadata_table.ExperimentDate{i}), frequency_tests_data_table.expdate));
            valid_threshold_rows = threshold_rows(~isnan([frequency_tests_data_table.thresholds_amplitude_mA(threshold_rows)]));
            freq_i = min([frequency_tests_data_table.frequencies_kHz(valid_threshold_rows)]);
            
            if strcmp(metadata_table.CuffType{i},'Monopolar')
                min_block_thresholds_monopolar = vertcat(min_block_thresholds_monopolar,freq_i);
            elseif strcmp(metadata_table.CuffType{i},'Bipolar')
                min_block_thresholds_bipolar = vertcat(min_block_thresholds_bipolar,freq_i);
            elseif strcmp(metadata_table.CuffType{i},'Tripolar')
                min_block_thresholds_tripolar = vertcat(min_block_thresholds_tripolar,freq_i);
            else
                error('Invalid metadata_table.CuffType value')
            end
        end
        
    end
end
fig_handle = compare_group_fitness({min_block_thresholds_monopolar,min_block_thresholds_bipolar,min_block_thresholds_tripolar});

set(fig_handle,'position',[435,606,233,233])
ylabel('min block frequency (kHz)')
xlabel('# contacts')
ylim([0 16])

