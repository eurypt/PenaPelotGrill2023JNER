%{

Plots figures 7,8, and 9 quantifying the effect of composite signals on
block threshold rms and and onset response.

%}

function plot_shape_effects_on_onset_auc()

cuff_type_to_plot = 'monopolar';
analyze(cuff_type_to_plot)

cuff_type_to_plot = 'bipolar';
analyze(cuff_type_to_plot)

cuff_type_to_plot = 'tripolar';
analyze(cuff_type_to_plot)


function analyze(cuff_type_to_plot)

% add path to internal dependencies (i.e., functions I wrote to facilitate
% plotting/analysis/etc.)
addpath('internal_dependencies');

%%% Load key data table
load('composite_signal_tests_data_table.mat','composite_signal_tests_data_table');

%%% Remove rows in which there was no valid block threshold
composite_signal_tests_data_table(isnan(composite_signal_tests_data_table.threshold_mA_RMS),:) = [];

%% Section 2: For each exp, make non-dome plots of the main metric (onset response or block threhsold) vs. inclination, and/or azimuth

% Specify plotting parameters and flags
%%% Specify vector of amplitudes at which to interpolate onset response;
%%% NOTE: if units_of_KHF_amplitude = 'normalized', then
%%% KHF_amplitude_to_plot specifies the value relative to block threshold
% if (~exist('amplitude_relative_to_BT','var') || isempty(amplitude_relative_to_BT))
amplitude_relative_to_BT = 1.1; %[1]; % [units specified by the next variable below]
% end
%sum(1.1>cellfun(@(x,y) max(x/y),composite_signal_tests_data_table_exp_i.KHF_amplitude_scatter_data,num2cell(composite_signal_tests_data_table_exp_i.threshold_mA)))

%%% Set to 'normalized' to normalize KHF by the block threshold amplitude;
%%% otherwise, set to 'mA'
% if (~exist('units_of_KHF_amplitude','var') || isempty(units_of_KHF_amplitude))
units_of_KHF_amplitude = 'normalized';
% end

%%% Set to 'onset' to plot onset response; set to 'threshold' to plot
%%% threshold; set to 'threshold_rms' to plot threshold in terms of RMS
main_metric_to_plot = 'threshold_rms';

%%% Set to 1 to plot difference between pairs (only applicable to datasets
%%% in which each waveform has a valid pair (i.e., pair_status = 'paired')
% if (~exist('flag_plot_pairwise_differences','var') || isempty(flag_plot_pairwise_differences))
flag_plot_pairwise_differences = 0;
% end



% Check plotting parameters
assert(any(strcmp(units_of_KHF_amplitude,{'normalized','mA'})),'Invalid value for units_of_KHF_amplitude');
assert(any(strcmp(main_metric_to_plot,{'threshold','threshold_rms','onset'})),'Invalid value for main_metric_to_plot');



% Identify unique datasets to process
if (flag_plot_pairwise_differences)
    is_pairwise = strcmp(composite_signal_tests_data_table.pair_status,'paired');
else
    is_pairwise = true(size(composite_signal_tests_data_table,1),1);
    unique_valid_experiment_dates = unique(composite_signal_tests_data_table.experiment_dates_all);
end

matching_cuff_type = strcmp(composite_signal_tests_data_table.cuff_type,cuff_type_to_plot);

unique_valid_experiment_dates = unique(composite_signal_tests_data_table.experiment_dates_all(is_pairwise & matching_cuff_type));



%- Store main_metric_to_plot (at the specified KHF amplitude), block
%- threshold, inclination, azimuth, and coefficients of each waveform
onset_response_cell_array = cell(length(unique_valid_experiment_dates),length(amplitude_relative_to_BT));
block_threshold_cell_array = cell(length(unique_valid_experiment_dates),length(amplitude_relative_to_BT));
inclination_cell_array = cell(length(unique_valid_experiment_dates),length(amplitude_relative_to_BT));
azimuth_cell_array = cell(length(unique_valid_experiment_dates),length(amplitude_relative_to_BT));
coefficients_cell_array = cell(length(unique_valid_experiment_dates),length(amplitude_relative_to_BT));
for i = 1:length(amplitude_relative_to_BT)
    
    %- For each unique experiment date...
    for exp_date_idx = 1:length(unique_valid_experiment_dates)
        
        %- Get the subset of data table relevant for that experiment
        exp_date_i = unique_valid_experiment_dates{exp_date_idx};
        composite_signal_tests_data_table_exp_i = composite_signal_tests_data_table(...
            strcmp(composite_signal_tests_data_table.experiment_dates_all,exp_date_i),:);
        
        %- For each test within this experiment...
        for row_idx = 1:size(composite_signal_tests_data_table_exp_i,1)
            
            %- Get each KHF amplitude within the present test; if specified,
            %- normalize KHF amplitude relative to block threshold
            if (strcmp(units_of_KHF_amplitude,'normalized'))
                KHF_amp_values_actual = composite_signal_tests_data_table_exp_i.KHF_amplitude_scatter_data{row_idx}/composite_signal_tests_data_table_exp_i.threshold_mA(row_idx);
            elseif (strcmp(units_of_KHF_amplitude,'mA'))
                KHF_amp_values_actual = composite_signal_tests_data_table_exp_i.KHF_amplitude_scatter_data{row_idx};
            else
                error('Invalid value of units_of_KHF_amplitude')
            end
            
            %- For each KHF amplitude within the present test, get the onset
            %- response
            normalized_AUC = composite_signal_tests_data_table_exp_i.AUC_onset_response_scatter_data{row_idx};
            
            %- Interpolate onset response at a fixed amplitude; remove
            %- non-unique KHF amplitudes for interpolation
            normalized_AUC_unique_KHF_amp = [];
            unique_KHF_amp = unique(KHF_amp_values_actual);
            for unique_amp_idx = 1:length(unique_KHF_amp)
                normalized_AUC_unique_KHF_amp(unique_amp_idx) = mean(normalized_AUC(KHF_amp_values_actual==unique_KHF_amp(unique_amp_idx)));
            end
            interpolated_AUC = interp1(unique_KHF_amp,normalized_AUC_unique_KHF_amp,amplitude_relative_to_BT(i),'linear',NaN);
            
            % If desired amplitude is greater than any of the sampled
            % amplitudes, extrapolate using nearest neighbors
            % extrapolation
            if (max(normalized_AUC_unique_KHF_amp)<amplitude_relative_to_BT(i))
                [~,max_idx] = max(normalized_AUC_unique_KHF_amp);
                interpolated_AUC = normalized_AUC_unique_KHF_amp(max_idx);
            end
            
            %- Calculate inclination and azimuth of waveform shape (i.e.,
            %- as characterized in f0 vs. 2f0 (sine) vs. 2f0 (cosine)
            %- space
            x = composite_signal_tests_data_table_exp_i.coefficients{row_idx}(2);
            y = composite_signal_tests_data_table_exp_i.coefficients{row_idx}(3);
            z = composite_signal_tests_data_table_exp_i.coefficients{row_idx}(1);
            theta = atan2d(sqrt(sum([x,y].^2,2)),z); % inclination
            phi = atan2d(y,x); % azimuth
            
            %- Store values
            onset_response_cell_array{exp_date_idx,i}(row_idx) = interpolated_AUC;
            block_threshold_cell_array{exp_date_idx,i}(row_idx) = composite_signal_tests_data_table_exp_i.threshold_mA(row_idx);
            inclination_cell_array{exp_date_idx,i}(row_idx) = theta;
            coefficients_cell_array{exp_date_idx,i}(row_idx,:) = composite_signal_tests_data_table_exp_i.coefficients{row_idx};
            
            %- If this dataset is f0 or 2*f0, set azimuth value to NaN,
            %- since there is no valid azimuth for non-mixed sinusoids
            ARBITRARY_PRECISION = 9;
            if (isequal(round(composite_signal_tests_data_table_exp_i.coefficients{row_idx},ARBITRARY_PRECISION),[1 0 0]))
                azimuth_cell_array{exp_date_idx,i}(row_idx) = NaN;
                %             plot(row_idx*[1 1],[0,max(cellfun(@max,composite_signal_tests_data_table_exp_i.KHF_amplitude_scatter_data))],'k--');
                %             hold on;
            elseif (isequal(round(composite_signal_tests_data_table_exp_i.coefficients{row_idx},ARBITRARY_PRECISION),[0 1 0]))
                azimuth_cell_array{exp_date_idx,i}(row_idx) = NaN;
                %             plot(row_idx*[1 1],[0,max(cellfun(@max,composite_signal_tests_data_table_exp_i.KHF_amplitude_scatter_data))],'r--');
                %             hold on;
            else
                azimuth_cell_array{exp_date_idx,i}(row_idx) = phi;
            end
            
            
            %{
        % if full or partial block occurred, plot onset response
        scatter3(thresholds *ones(size(KHF_amp_values_actual(full_or_partial_block_occurred))),...
            KHF_amp_values_actual(full_or_partial_block_occurred),normalized_AUC(full_or_partial_block_occurred),...
            50,normalized_AUC(full_or_partial_block_occurred),'fill','MarkerEdgeColor','k')
        hold on
        
        % otherwise, just plot a hollow scatter point
        scatter3(thresholds *ones(size(KHF_amp_values_actual(~full_or_partial_block_occurred))),...
            KHF_amp_values_actual(~full_or_partial_block_occurred),normalized_AUC(~full_or_partial_block_occurred),...
            50,'MarkerEdgeColor','k')
%         caxis([0 CUTOFF_FOR_PLOTTING])
        % ylim([0 max(KHF_amp_values_actual)])
        hold on;
            %}
        end
    end
end


%- Specify main metric to plot
if (strcmp(main_metric_to_plot,'onset'))
    main_metric_cell_array = onset_response_cell_array;
    side_metric_cell_array = block_threshold_cell_array;
    main_label_str = 'onset response (N*s)';
    side_label_str = 'block threshold (mA)';
elseif (strcmp(main_metric_to_plot,'threshold'))
    main_metric_cell_array = block_threshold_cell_array;
    side_metric_cell_array = onset_response_cell_array;
    main_label_str = 'block threshold (mA)';
    side_label_str = 'onset response (N*s)';
elseif (strcmp(main_metric_to_plot,'threshold_rms'))
    main_metric_cell_array = block_threshold_cell_array;
    for i = 1:length(main_metric_cell_array)
        for j = 1:length(main_metric_cell_array{i})
            main_metric_cell_array{i}(j) = main_metric_cell_array{i}(j)*rms(construct_waveform(10e3,coefficients_cell_array{i}(j,:)));
        end
    end
    
    side_metric_cell_array = onset_response_cell_array;
    main_label_str = 'block threshold (mA RMS)';
    side_label_str = 'onset response (N*s)';
else
    error('invalid value for main_metric_to_plot')
end


%% Section 3: Rank waveforms by onset response and threshold (equally weighted) and plot (aligned)

if (~strcmp(main_metric_to_plot,'onset')) % skip the waveforms plotting if metric being plotted is onset response
    
    %%% Specify whether to sort with 'ascend' (i.e., first few waveforms
    %%% will be the best) or with 'descend' (i.e., first few waveforms will
    %%% be the worst)
    sort_type_all = {'ascend','descend'};
    
    % For each dataset, reconstruct the waveforms, sort them, and plot the top
    % N waveforms (aligned) in a single plot
        all_top_waveforms = [];
    for sort_type_idx = 1:length(sort_type_all)
        sort_type = sort_type_all{sort_type_idx};
        
        for i = 1:length(unique_valid_experiment_dates)
            if(~strcmp(main_metric_to_plot,'threshold_rms'))
                warning('Waveform efficiency plots should not be performed for the ''threshold'' metric; use ''threshold_rms'' only...'); pause(3)
            end
            
            % Reconstruct waveforms in time domain and align the top N waveforms
            waveforms_i = [];
            coefficients_all = coefficients_cell_array{i};
            % only include sum of sine waveforms
            ARBITRARY_PRECISION = 9;
            all_evaluated_inputs_i = coefficients_all;
            base_freq_indices = find(round(all_evaluated_inputs_i(:,1),ARBITRARY_PRECISION)==1);
            twice_base_freq_indices = find(round(all_evaluated_inputs_i(:,2),ARBITRARY_PRECISION)==1 | ...
                round(all_evaluated_inputs_i(:,3),ARBITRARY_PRECISION)==1);
            other_input_indices = setdiff(1:length(all_evaluated_inputs_i),[base_freq_indices; twice_base_freq_indices]);
            coefficients_subset = coefficients_all(other_input_indices,:);
            n_waveforms_i = size(coefficients_subset,1);
            
            for j = 1:n_waveforms_i
                % Reconstruct the waveform at some dummy (irrelevant) frequency
                dummy_freq = 10e3; % [Hz]
                [custom_waveform_ij,~] = construct_waveform(dummy_freq,coefficients_subset(j,:));
                waveforms_i = vertcat(waveforms_i,custom_waveform_ij);
            end
            
            % Ignore waveforms that had no valid block threshold or no valid
            % onset response value
            block_thresholds_i = main_metric_cell_array{i}(other_input_indices);
            onset_responses_i = onset_response_cell_array{i}(other_input_indices);
            indices_to_ignore = find(isnan(block_thresholds_i) | isnan(onset_responses_i));
            waveforms_i(indices_to_ignore,:) = [];
            block_thresholds_i(indices_to_ignore) = [];
            onset_responses_i(indices_to_ignore) = [];
            
            % Sort waveforms by onset response and threshold (equally weighted) (either
            % ascending or descending)
            %%% Specify weighting (unitless, from 0 to 1); 0 will use block
            %%% threshold only; 1 will use onset response only
            weighting = 0;
            %%% Sort waveforms
            multi_objective_values = (1-weighting)*block_thresholds_i + weighting*onset_responses_i;
            [sorted_multi_objective_values,sort_order] = sort(multi_objective_values,sort_type);
            sorted_waveforms_i = waveforms_i(sort_order,:);
            
            
            % Align the top N waveforms to the first waveform, where N is user-specified
            %%% Specify align method
            align_method = 2;
            %%%
            N = 7;
            top_N_waveforms = sorted_waveforms_i(1:N,:);
            all_top_waveforms(i,1:size(top_N_waveforms,2),sort_type_idx) = top_N_waveforms(1,:);
            
        end
        
    end
    
    % Align the best waveforms across all nerves and the inverse of the
    % worst waveforms across all nerves; flip the sign for alinging then
    % unflip it for storing back in the original variable
    temp = alignsignals_circ([all_top_waveforms(:,:,1); -all_top_waveforms(:,:,2)]);
    all_top_waveforms(:,:,1) = temp(1:length(unique_valid_experiment_dates),:);
    all_top_waveforms(:,:,2) = -temp(length(unique_valid_experiment_dates) + (1:length(unique_valid_experiment_dates)),:);
    
    % Plot all the best waveforms across all nerves and (in a separate
    % window) all the worst waveforms across all nerves
    for sort_type_idx = 1:length(sort_type_all)
        
        figure('color',[1 1 1]); plot(all_top_waveforms(:,:,sort_type_idx)','LineWidth',2);
        hold on;
        
        set(gca,'XTickLabel','','FontSize',14);
        ylim([-1 1])
        set(gca,'YTick','');
        
        % Add a zero line
        xlim_i = get(gca,'XLim');
        plot(xlim_i,[0 0],'-','color',[0.5 0.5 0.5])
        
        legend(unique_valid_experiment_dates)
    end
    
    % Quantify and print the zero-crossing width of the largest positive
    % positive or negative peak

    for sort_type_idx = 1:length(sort_type_all)
        % Find the max abs peak and the peak sign
        [~,max_idx] = max(abs(all_top_waveforms(:,:,sort_type_idx)),[],2);
        peak_sign = zeros(size(max_idx));
        for j = 1:size(all_top_waveforms,1)
            peak_sign(j) = sign(all_top_waveforms(j,max_idx(j),sort_type_idx));
        end
        
        % Find the size of the peak and the oppositely singed peak
        peaks = [];
        for j = 1:size(all_top_waveforms,1)
            peaks(j,1) = max(peak_sign(j)*all_top_waveforms(j,:,sort_type_idx),[],2);
            peaks(j,2) = max(-1*peak_sign(j)*all_top_waveforms(j,:,sort_type_idx),[],2);
        end
        
        % Replicate the signal twice; adjust the peak location to capture the
        % middle wavefrom's peak
        replicated_signal = repmat(all_top_waveforms(:,:,sort_type_idx),1,3);
        
        crossings = [];
        
        % Find the last zero cross before the peak
        for j = 1:size(all_top_waveforms,1)
            peak_loc = max_idx(j)+size(all_top_waveforms,2);
            crossings(j,1) = find(peak_sign(j)*replicated_signal(j,1:peak_loc)<0,1,'last');
        end
        
        % Find the first zero crossing after the peak
        for j = 1:size(all_top_waveforms,1)
            peak_loc = max_idx(j)+size(all_top_waveforms,2);
            crossings(j,2) = find(peak_sign(j)*replicated_signal(j,peak_loc:end)<0,1,'first')+(peak_loc-1);
        end
        
        % Calculate the zero-crossing duration
        zero_cross_duration_percent = 100*diff(crossings,[],2)/size(all_top_waveforms,2); % [percent]
        
        % Print for all biggest peak zero-crossing durations
        fprintf('percent zero cross duration relative to period\n');
        fprintf('%d\n',round(zero_cross_duration_percent,1))
        
        
        % Print for all peak ratios
        peark_ratios = peaks(:,1)./peaks(:,2);
        fprintf('peak_ratios (larger to smaller)\n');
        fprintf('%0.3f\n',peark_ratios)
    end
end

%% Section 4: undescribed for now
%- If specified to plot pairwise data, replace the data with pairwise
%- differences
if (flag_plot_pairwise_differences)
    % for each element, if it has a valid azimuth value, find its pair, then
    % store the difference between the onset response pair into a 'paired onset
    % response' array, but only for the larger of the azimuth values; for the
    % smaller of the azimuth values, store NaN into the 'paired onset response'
    % array
    paired_main_metric_cell_array = cell(size(main_metric_cell_array));
    for i = 1:size(azimuth_cell_array,1)
        for j = 1:size(azimuth_cell_array,2)
            for k = 1:length(azimuth_cell_array{i,j})
                if ~isnan(azimuth_cell_array{i,j}(k)) % is there is a valid azimuth value (since pure sinusoids do *not* have a valid azimuth value)
                    
                    % find the pair
                    ARBITRARY_PRECISION = 6;
                    k_pair_idx = find(...
                        (round(azimuth_cell_array{i,j},ARBITRARY_PRECISION)==round(azimuth_cell_array{i,j}(k)-180,ARBITRARY_PRECISION) | ...
                        round(azimuth_cell_array{i,j},ARBITRARY_PRECISION)==round(azimuth_cell_array{i,j}(k)+180,ARBITRARY_PRECISION)) &...
                        round(inclination_cell_array{i,j},ARBITRARY_PRECISION)==round(inclination_cell_array{i,j}(k),ARBITRARY_PRECISION));
                    assert(length(k_pair_idx)==1,'one and only one pair must be found for each test...')
                    
                    % store the difference if the current element's azimuth is
                    % bigger than its pair's azimuth
                    if (azimuth_cell_array{i,j}(k)>azimuth_cell_array{i,j}(k_pair_idx))
                        paired_main_metric_cell_array{i,j}(k) = main_metric_cell_array{i,j}(k)-main_metric_cell_array{i,j}(k_pair_idx);
                    else
                        paired_main_metric_cell_array{i,j}(k) = NaN;
                    end
                    
                end
            end
        end
    end
    
    main_metric_cell_array = paired_main_metric_cell_array;
    main_label_str = {'pairwise difference in',main_label_str};
end

%{
%- Initialize figure window
% window_width_factor = (125/15); % plot window width will be adaptively expanded to window_width_factor times the number of tests to be plotted
% window_width = window_width_factor*size(composite_signal_tests_data_table_exp_i,1);
%}

%- Plot all combinations of variables against each other
plot_combos = 3;
for i = 1:length(amplitude_relative_to_BT)
    
    % Set up label for onset response if applicable
    if strcmp(main_metric_to_plot,'onset')
        amp_suffix_str = sprintf('%0.2f*BT',amplitude_relative_to_BT(i));
    else
        amp_suffix_str = '';
    end
    
    for combo_idx = 1:plot_combos
        figure('color',[1 1 1],'position',[42,171,1370,222]);
        for exp_date_idx = 1:length(unique_valid_experiment_dates)
            
            
            
            %- Make the plot of the specified combo
            switch combo_idx
                case 1 % Plot onset response vs. block threshold
                    %- Assign subplot panel
                    subplot(1,length(unique_valid_experiment_dates),exp_date_idx);
                    
                    %- Plot onset response at a fixed KHF amplitude vs. block
                    %- threshold
                    scatter(side_metric_cell_array{exp_date_idx,i},main_metric_cell_array{exp_date_idx,i},25,'k');
                    ylabel(main_label_str)
                    xlabel(side_label_str)
                    title({[unique_valid_experiment_dates{exp_date_idx},' ',sprintf('n=%d',sum(~isnan(main_metric_cell_array{exp_date_idx,i})))],amp_suffix_str})
                    set(gca,'FontSize',10)
                    
                case 2 % Plot onset response at a fixed KHF amplitude vs. inclination
                    if (flag_plot_pairwise_differences)
                        subplot(1,length(unique_valid_experiment_dates),exp_date_idx);
                    else
                        h = subplot(1,3*length(unique_valid_experiment_dates),(exp_date_idx-1)*3+[1 2]);
                    end
                    
                    % sine sum
                    
                    % plot onset vs. inclination, with azimuth as color
                    scatter(inclination_cell_array{exp_date_idx,i},...
                        main_metric_cell_array{exp_date_idx,i},25,...
                        'k','LineWidth',1);
                    
                    %{
                        % plot onset vs. inclination, with azimuth as color
                        scatter(inclination_cell_array{exp_date_idx,i},...
                        onset_response_cell_array{exp_date_idx,i},25,...
                        azimuth_cell_array{exp_date_idx,i},'LineWidth',3);
                        %}
                        xlabel('inclination (deg)')
                        if (exp_date_idx==1) % only label the leftmost panel
                            ylabel(main_label_str)
                        end
                        set(gca,'FontSize',10)
                        title({[unique_valid_experiment_dates{exp_date_idx},' ',sprintf('n=%d',sum(~isnan(main_metric_cell_array{exp_date_idx,i})))],amp_suffix_str})
                        
                        
                        if (~flag_plot_pairwise_differences)
                            is_f0_or_2f0 = isnan(azimuth_cell_array{exp_date_idx,i});
                            %- Plot onset response at a fixed KHF amplitude vs.
                            %- azimuth; set inclination as the color
                            subplot(1,3*length(unique_valid_experiment_dates),(exp_date_idx-1)*3+3);
                            % pure sine
                            
                            scatter(inclination_cell_array{exp_date_idx,i}(is_f0_or_2f0),...
                                main_metric_cell_array{exp_date_idx,i}(is_f0_or_2f0),25,...
                                'k','filled');
                            ylim(get(h,'YLim'))
                            set(gca,'YTickLabel','')
                            base_frequency_Hz = unique(composite_signal_tests_data_table.base_freq_Hz(strcmp(unique_valid_experiment_dates{exp_date_idx},composite_signal_tests_data_table.experiment_dates_all)));
%                             set(gca,'XTick',[0 90],'XTickLabel',{'f_0','2f_0'})
                            f0_label = sprintf('%d',round(base_frequency_Hz/1e3));
                            twice_f0_label = sprintf('%d',2*round(base_frequency_Hz/1e3));
                            set(gca,'XTick',[0 90],'XTickLabel',{f0_label,twice_f0_label})
                            xlim([-45 135]);
                            xlabel('(kHz)')
                            set(gca,'FontSize',10)
                        end
                        
                case 3 % Plot onset response vs. azimuth, using inclination as color
                    %- Identify the f0 and 2f0 waveforms so they can be
                    %- plotted separately from the sinesum
                    is_f0_or_2f0 = isnan(azimuth_cell_array{exp_date_idx,i});
                    f0_or_2f0_indices = find(is_f0_or_2f0);
                    composite_signal_indices = find(~is_f0_or_2f0);
                    
                    warning('ignoring the very first f0 trial...'); pause(3)
                    f0_or_2f0_indices(1) = [];
                    
                    if (~flag_plot_pairwise_differences)
                        
                        %- Plot onset response at a fixed KHF amplitude vs.
                        %- azimuth; set inclination as the color
                        
                        % fixed so that all the data is plotted with the
                        % same numebr of columns; this will make all the
                        % plots the same width
                        NUMBER_OF_PLOTS = 7; % length(unique_valid_experiment_dates
                        
                        g = subaxis(1,3*NUMBER_OF_PLOTS,(exp_date_idx-1)*3+3,'ML',0.05,'MT',0.2,'MB',0.2,'MR',0,'SpacingHoriz',0);
                        % pure sine
                        
                        % Plot the f0 and 2f0 data in pairs such that a
                        % line connects the pairs, since they were all
                        % sampled in pairs (in randomized order)
                        xdata = inclination_cell_array{exp_date_idx,i};
                        ydata = main_metric_cell_array{exp_date_idx,i};
                        idx = 1;
                        while idx <= length(f0_or_2f0_indices)
%                         for pair_idx = 1:length(xdata)
                            if (idx < length(f0_or_2f0_indices) && f0_or_2f0_indices(idx+1)==(f0_or_2f0_indices(idx)+1))
                                plot(xdata(f0_or_2f0_indices(idx) + [0,1]),...
                                    ydata(f0_or_2f0_indices(idx) + [0,1]), 'k.-', 'MarkerSize', 20);
                                idx = idx + 2;
                            else
                                plot(xdata(f0_or_2f0_indices(idx)),...
                                    ydata(f0_or_2f0_indices(idx)), 'k.-', 'MarkerSize', 20);
                                idx = idx + 1;
                            end
                            hold on;
                        end
                        box off
                        
%                         scatter(inclination_cell_array{exp_date_idx,i}(f0_or_2f0_indices),...
%                             main_metric_cell_array{exp_date_idx,i}(f0_or_2f0_indices),25,...
%                             'k','filled');
                        %                         ylim(get(h,'YLim'))
                        ylim([0 1.1*max(main_metric_cell_array{exp_date_idx,i})])
                        set(gca,'YTickLabel','')
                        base_frequency_Hz = unique(composite_signal_tests_data_table.base_freq_Hz(strcmp(unique_valid_experiment_dates{exp_date_idx},composite_signal_tests_data_table.experiment_dates_all)));
                        %                             set(gca,'XTick',[0 90],'XTickLabel',{'f_0','2f_0'})
                        f0_label = sprintf('%d',round(base_frequency_Hz/1e3));
                        twice_f0_label = sprintf('%d',2*round(base_frequency_Hz/1e3));
                        set(gca,'XTick',[0 90],'XTickLabel',{f0_label,twice_f0_label})
                        xlim([-45 135]);
                        xlabel('(kHz)')
                        set(gca,'FontSize',10)
                        
                        
                    end
                    
                    %- Plot onset response at a fixed KHF amplitude vs.
                    %- azimuth; set inclination as the color
                    if (flag_plot_pairwise_differences)
                        subplot(1,length(unique_valid_experiment_dates),exp_date_idx);
                    else
                        h = subaxis(1,3*NUMBER_OF_PLOTS,(exp_date_idx-1)*3+[1 2]); %,'ML',0,'MT',0,'MB',0,'MR',0,'SpacingHoriz',0);
                    end
                    
                    
                    % sine sum (plotting as scatter3 mainly so that the
                    % data cursor displays the inclination too, which makes
                    % it easier to find the original raw datasets)
                    scatter3(azimuth_cell_array{exp_date_idx,i}(composite_signal_indices),main_metric_cell_array{exp_date_idx,i}(~is_f0_or_2f0),inclination_cell_array{exp_date_idx,i}(~is_f0_or_2f0),...
                        35,inclination_cell_array{exp_date_idx,i}(composite_signal_indices),'filled'); %,'LineWidth',1);
                    view([0 90]); % set the default view as XY
                    
                    xlim(215*[-1 1]);
                    set(gca,'XTick',-180:180:180)
                    
                    xlabel('azimuth (deg)')
                    if (exp_date_idx==1) % only label the leftmost panel
                        ylabel(main_label_str)
                    end
                    set(gca,'FontSize',10)
                    if (~flag_plot_pairwise_differences)
                        set(gca,'YLim',get(g,'YLim'))
                        linkaxes([g,h],'y')
                    end
                    title({[unique_valid_experiment_dates{exp_date_idx},' ',sprintf('n=%d',sum(~isnan(main_metric_cell_array{exp_date_idx,i})))],amp_suffix_str})
                    
                otherwise
                    error('Invalid combo_idx')
            end
            
        end
        
        flag_show_colorbar = 0;
        if (combo_idx==3 && flag_show_colorbar)
            h = colorbar;
            ylabel(h,'inclination (deg)')
        end
    end
end
    
    %% Section 9: Calculate the onset response
    
    %%% Create storage table for the all KHF amplitudes and onset responses
    %%% as well as the single-amplitude interpolated onset response; this
    %%% storage table will have all rows corresponding to those in
    %%% master_thresholds_table and will ultimately be joined with it
    onset_response_table = table(nan(size(composite_signal_tests_data_table.threshold_mA)),...
        cell(size(composite_signal_tests_data_table.threshold_mA)),...
        cell(size(composite_signal_tests_data_table.threshold_mA)),...
        cell(size(composite_signal_tests_data_table.threshold_mA)),...
        'VariableNames',{'onset_response_at_single_amp_Ns','onset_responses_all',...
        'KHF_amplitudes_all','normalized_KHF_amplitudes_all'});
    
    %  For every unique test for which onset response is quantified, iterate to
    %  get the AUC vs. KHF amp curve. If that experiment with frequency combo
    %  has a block threshold available, normalize the curve to the block
    %  threshold, and store the AUC vs. raw amplitudes curve as well as the AUC
    %  vs. normalized amplitudes curve. Then, store as 'threshold' the AUC at a
    %  fixed percentage of block threshold
    % unique_exp = unique(composite_signal_tests_data_table.experiment_dates_all);
    % for exp_idx = 1:length(unique_exp)
    %     exp_i = unique_exp{exp_idx};
    %     test_row_indices = find(strcmp(composite_signal_tests_data_table.experiment_dates_all,exp_i));
    %     for test_row_idx = 1:length(test_row_indices)
    for match_idx = 1:size(composite_signal_tests_data_table,1)
        %         row_ij = test_row_indices(test_row_idx);
        
        AUC_ij = composite_signal_tests_data_table.AUC_onset_response_scatter_data{match_idx};
        KHF_ij = composite_signal_tests_data_table.KHF_amplitude_scatter_data{match_idx};
        
        % combine repeated amplitude points as mean
        ARBITRARY_PRECISION = 8;
        [~,~,IA] = unique(round(KHF_ij,ARBITRARY_PRECISION));
        KHF_ij_merge_repeats = zeros(max(IA),1);
        AUC_ij_merge_repeats = zeros(max(IA),1);
        for k = 1:max(IA)
            KHF_ij_merge_repeats(k) = mean(KHF_ij(IA==k));
            AUC_ij_merge_repeats(k) = mean(AUC_ij(IA==k));
        end
        
        % check whether this experiment has a block threshold available
        %         coefficients_all = cell2mat(composite_signal_tests_data_table.coefficients);
        %         match_idx = find(strcmp(composite_signal_tests_data_table.experiment_dates_all{row_ij},composite_signal_tests_data_table.experiment_dates_all) & ...
        %             composite_signal_tests_data_table.coefficients{row_ij}(1)==coefficients_all(:,1) & ...
        %             composite_signal_tests_data_table.coefficients{row_ij}(2)==coefficients_all(:,2) & ...
        %             composite_signal_tests_data_table.coefficients{row_ij}(3)==coefficients_all(:,3));
        if (~isempty(match_idx))
            
            % Get the block threshold
            BT_ij = composite_signal_tests_data_table.threshold_mA(match_idx);
            
            % Normalize KHF amplitudes by BT
            normalized_KHF_ij_merge_repeats = KHF_ij_merge_repeats/BT_ij;
            
            % Interpolate the onset response at the desired amplitude
            % (i.e., at amplitude_relative_to_BT)
            onset_response_table.onset_response_at_single_amp_Ns(match_idx) = interp1(normalized_KHF_ij_merge_repeats,AUC_ij_merge_repeats,amplitude_relative_to_BT);
            
            % If desired amplitude is greater than any of the sampled
            % amplitudes, extrapolate using nearest neighbors
            % extrapolation
            if (max(normalized_KHF_ij_merge_repeats)<amplitude_relative_to_BT)
                [~,max_idx] = max(normalized_KHF_ij_merge_repeats);
                onset_response_table.onset_response_at_single_amp_Ns(match_idx) = AUC_ij_merge_repeats(max_idx);
            end
            
            % If interpolated onset response value is negative, set it to positive
            % and print out a warning plus the value of the
            if (onset_response_table.onset_response_at_single_amp_Ns(match_idx)<0)
                warning('setting the following onset response value to be positive: %0.8f',onset_response_table.onset_response_at_single_amp_Ns(match_idx))
                onset_response_table.onset_response_at_single_amp_Ns(match_idx) = abs(onset_response_table.onset_response_at_single_amp_Ns(match_idx));
                
                % pause briefly so user can see the warning
                pause(2);
            end
            
            % Store the AUC vs. raw amplitudes curve as well as the
            % AUC vs. normalized amplitudes curve
            onset_response_table.onset_responses_all{match_idx,1} = AUC_ij_merge_repeats;
            onset_response_table.KHF_amplitudes_all{match_idx,1} = KHF_ij_merge_repeats;
            onset_response_table.normalized_KHF_amplitudes_all{match_idx,1} = normalized_KHF_ij_merge_repeats;
        end
        %     end
    end
    
    %%% Add the onset response at single amplitude, AUC, raw amplidues, and
    %%% normalized amplitudes to the thresholds table
    master_thresholds_table = cat(2,composite_signal_tests_data_table,onset_response_table);
    
    %%% Again, remove rows that do not have a valid threshold value
    assert(~any(isnan(master_thresholds_table.onset_response_at_single_amp_Ns)));
    %{
if (strcmp(metric_str,'onset'))
    composite_signal_tests_data_table(isnan(composite_signal_tests_data_table.onset_response_at_single_amp_Ns),:) = [];
end
    %}
    
    %%% Plot the max amplitudes sampled relative to BT for every data point
    figure('color',[1 1 1]); histogram(cellfun(@(x) max(x),master_thresholds_table.normalized_KHF_amplitudes_all))
    xlabel({'scaling factor value','(i.e., amplitude relative to BT)'})
    ylabel('count')
    set(gca,'FontSize',14)
    
    
    %% Section 10: Compare (via scatter plots containing points from all datasets) most efficient vs. least efficient vs. f0 vs 2f0 in terms of BT and onset response
    
    % Specify whether to use the mean, median, min, or max of replicate values
    
    % % Select the metric that will be used for subsequent plotting and
    % % analysis
    % metric_str = 'power';
    % metric = zeros(size(master_thresholds_table,1),1);
    % if (strcmp(metric_str,'current'))
    %     metric = master_thresholds_table.threshold_mA;
    %     thresholds_str = 'threshold amplitude (mA)';
    % elseif (strcmp(metric_str,'power'))
    %     metric = master_thresholds_table.threshold_mA_RMS;
    %     thresholds_str = 'threshold amplitude (mW)';
    % elseif (strcmp(metric_str,'onset'))
    %     metric = master_thresholds_table.onset_response_at_single_amp_Ns;
    %     thresholds_str = 'onset response (N*s)';
    % end
    
    % Get the values
    
    % minimum_frequencies_for_each_nerve = cellfun(@(x) min(...
    %     master_thresholds_table.frequencies_kHz(strcmp(master_thresholds_table.experiment_dates_all,x))), unique(master_thresholds_table.experiment_dates_all));
    
    % coefficients_all = master_thresholds_table.coefficients
    % energy_at_minimum_frequences_for_each_nerve = cellfun(@(x,y) metric(...
    %     ==x & strcmp(master_thresholds_table.experiment_dates_all,y)), ...
    %     num2cell(minimum_frequencies_for_each_nerve), unique(master_thresholds_table.experiment_dates_all));
    
    % energy_at_minimum_energy_for_each_nerve = cellfun(@(x) min(...
    %     metric(strcmp(master_thresholds_table.experiment_dates_all,x))), unique(master_thresholds_table.experiment_dates_all));
    
    % Identify unique nerves
    unique_nerves = unique(master_thresholds_table.experiment_dates_all);
    
    % For each nerve, find the minimum frequency and its energy and onset
    % response
    min_freq_data_table = table(cell(length(unique_nerves),1),zeros(length(unique_nerves),1),zeros(length(unique_nerves),1),...
        'VariableNames',{'coefficients','threshold_mA_RMS','onset_response_Ns'});
    for i = 1:length(unique_nerves)
        test_row_indices = find(strcmp(master_thresholds_table.experiment_dates_all,unique_nerves{i}));
        coefficient_min_freq = [1 0 0];
        min_freq_data_table.coefficients{i} = coefficient_min_freq;
        matching_rows = find(cellfun(@(x) x(1)==coefficient_min_freq(1) && x(2)==coefficient_min_freq(2) && x(3)==coefficient_min_freq(3), master_thresholds_table.coefficients(test_row_indices)));
        min_freq_data_table.threshold_mA_RMS(i) = median(master_thresholds_table.threshold_mA_RMS(test_row_indices(matching_rows)));
        min_freq_data_table.onset_response_Ns(i) = median(master_thresholds_table.onset_response_at_single_amp_Ns(test_row_indices(matching_rows)));
    end
    
    
    % For each nerve, find the maximum frequency and its energy and onset
    % response
    max_freq_data_table = table(cell(length(unique_nerves),1),zeros(length(unique_nerves),1),zeros(length(unique_nerves),1),...
        'VariableNames',{'coefficients','threshold_mA_RMS','onset_response_Ns'});
    for i = 1:length(unique_nerves)
        test_row_indices = find(strcmp(master_thresholds_table.experiment_dates_all,unique_nerves{i}));
        coefficient_max_freq = [0 1 0];
        max_freq_data_table.coefficients{i} = coefficient_max_freq;
        matching_rows = find(cellfun(@(x) x(1)==coefficient_max_freq(1) && x(2)==coefficient_max_freq(2) && x(3)==coefficient_max_freq(3), master_thresholds_table.coefficients(test_row_indices)));
        max_freq_data_table.threshold_mA_RMS(i) = median(master_thresholds_table.threshold_mA_RMS(test_row_indices(matching_rows)));
        max_freq_data_table.onset_response_Ns(i) = median(master_thresholds_table.onset_response_at_single_amp_Ns(test_row_indices(matching_rows)));
    end
    
    
    % For each nerve, find the minimum energy and the frequency and onset
    % response at the minimum energy
    min_energy_data_table = table(cell(length(unique_nerves),1),zeros(length(unique_nerves),1),zeros(length(unique_nerves),1),...
        'VariableNames',{'coefficients','threshold_mA_RMS','onset_response_Ns'});
    for i = 1:length(unique_nerves)
        test_row_indices = find(strcmp(master_thresholds_table.experiment_dates_all,unique_nerves{i}) & ...
            ~(cellfun(@(x) x(1)==coefficient_min_freq(1) && x(2)==coefficient_min_freq(2) && x(3)==coefficient_min_freq(3), master_thresholds_table.coefficients)) & ...
            ~(cellfun(@(x) x(1)==coefficient_max_freq(1) && x(2)==coefficient_max_freq(2) && x(3)==coefficient_max_freq(3), master_thresholds_table.coefficients)));
        [min_energy_data_table.threshold_mA_RMS(i),idx] = min(master_thresholds_table.threshold_mA_RMS(test_row_indices));
        min_energy_data_table.coefficients{i} = master_thresholds_table.coefficients{test_row_indices(idx)};
        min_energy_data_table.onset_response_Ns(i) = master_thresholds_table.onset_response_at_single_amp_Ns(test_row_indices(idx));
    end
    
    
    % For each nerve, find the minimum energy and the frequency and onset
    % response at the minimum energy
    max_energy_data_table = table(unique_nerves,cell(length(unique_nerves),1),zeros(length(unique_nerves),1),zeros(length(unique_nerves),1),...
        'VariableNames',{'unique_nerve','coefficients','threshold_mA_RMS','onset_response_Ns'});
    for i = 1:length(unique_nerves)
        test_row_indices = find(strcmp(master_thresholds_table.experiment_dates_all,unique_nerves{i}) & ...
            ~(cellfun(@(x) x(1)==coefficient_min_freq(1) && x(2)==coefficient_min_freq(2) && x(3)==coefficient_min_freq(3), master_thresholds_table.coefficients)) & ...
            ~(cellfun(@(x) x(1)==coefficient_max_freq(1) && x(2)==coefficient_max_freq(2) && x(3)==coefficient_max_freq(3), master_thresholds_table.coefficients)));
        [max_energy_data_table.threshold_mA_RMS(i),idx] = max(master_thresholds_table.threshold_mA_RMS(test_row_indices));
        max_energy_data_table.coefficients{i} = master_thresholds_table.coefficients{test_row_indices(idx)};
        max_energy_data_table.onset_response_Ns(i) = master_thresholds_table.onset_response_at_single_amp_Ns(test_row_indices(idx));
    end
    
    
    % Get electrode identity to group the plots by electrode type
    electrode_types = cellfun(@(x) ...
        master_thresholds_table.cuff_type{find(strcmp(master_thresholds_table.experiment_dates_all,x))}, unique_nerves,...
        'UniformOutput',false);
    
    % Plot the various combos
    xy_combos = {
        max_freq_data_table.threshold_mA_RMS,min_freq_data_table.threshold_mA_RMS, 'thresh (mA_{RMS}) at 2f_0','thresh (mA_{RMS}) at f_0'
        min_freq_data_table.onset_response_Ns, max_freq_data_table.onset_response_Ns,'onset resp (N*sec) at f_0','onset resp (N*sec) at 2f_0'
        min_energy_data_table.threshold_mA_RMS,min_freq_data_table.threshold_mA_RMS, 'thresh (mA_{RMS}) at most efficient composite','thresh (mA_{RMS}) at f_0'
        min_freq_data_table.onset_response_Ns, min_energy_data_table.onset_response_Ns,'onset resp (N*sec) at f_0','onset resp (N*sec) at most efficient composite'
        max_energy_data_table.threshold_mA_RMS,min_freq_data_table.threshold_mA_RMS, 'thresh (mA_{RMS}) at least efficient composite','thresh (mA_{RMS}) at f_0'
        min_freq_data_table.onset_response_Ns, max_energy_data_table.onset_response_Ns,'onset resp (N*sec) at f_0','onset resp (N*sec) at least efficient composite'
        min_energy_data_table.threshold_mA_RMS,max_freq_data_table.threshold_mA_RMS, 'thresh (mA_{RMS}) at most efficient composite','thresh (mA_{RMS}) at 2f_0'
        max_freq_data_table.onset_response_Ns, min_energy_data_table.onset_response_Ns,'onset resp (N*sec) at 2f_0','onset resp (N*sec) at most efficient composite'
        max_energy_data_table.threshold_mA_RMS,max_freq_data_table.threshold_mA_RMS, 'thresh (mA_{RMS}) at least efficient composite','thresh (mA_{RMS}) at 2f_0'
        max_freq_data_table.onset_response_Ns, max_energy_data_table.onset_response_Ns,'onset resp (N*sec) at 2f_0','onset resp (N*sec) at least efficient composite'
        max_energy_data_table.threshold_mA_RMS,min_energy_data_table.threshold_mA_RMS,'thresh (mA_{RMS}) at least efficient composite','thresh (mA_{RMS}) at most efficient composite'
        min_energy_data_table.onset_response_Ns, max_energy_data_table.onset_response_Ns,'onset resp (N*sec) at most efficient composite','onset resp (N*sec) at least efficient composite'
        % max freq vs. min energy (thresholds)
        % max freq vs. min energy (onset response)
        %     min_freq_data_table.threshold_mA_RMS, min_onset_data_table.threshold_mA_RMS,'thresh at min freq','thresh at min onset resp'
        %     min_freq_data_table.onset_response_Ns, min_onset_data_table.onset_response_Ns,'onset resp at min freq','onset resp at min onset resp'
        %     min_freq_data_table.frequency_kHz, min_energy_data_table.frequency_kHz,'freq at min freq','freq at most efficient composite'
        %     min_freq_data_table.frequency_kHz, min_onset_data_table.frequency_kHz,'freq at min freq','freq at min onset resp'
        };
    
    for i = 1:size(xy_combos,1)
        x = xy_combos{i,1};
        y = xy_combos{i,2};
        xlabel_str = xy_combos{i,3};
        ylabel_str = xy_combos{i,4};
        
        figure_title_str = sprintf('%s (x) vs. %s (y)',xlabel_str,ylabel_str);
        custom_figure('Name',figure_title_str,'NumberTitle','off','position',[680,558,459,420]);
        group_labels = {'monopolar','bipolar','tripolar'};
        group_colors = {[     0    0.4470    0.7410],[0.8500    0.3250    0.0980],[0.9290    0.6940    0.1250]};
        group_data = cellfun(@(group_i_label) ...
            [x(strcmp(electrode_types,group_i_label)),y(strcmp(electrode_types,group_i_label))],...
            group_labels,'UniformOutput',false);
        scatter_pairwise_comparison(group_data,group_colors,xlabel_str,ylabel_str,0,group_labels)
        xlabel('')
        ylabel('')
%         legend({'monopolar','bipolar','tripolar'},'Location','NorthWest')
        
        % Perform stat tests
        p = signrank(x,y);
        fprintf('%s: p = %0.6e\n',figure_title_str,p);
        set(gca,'FontSize',20)
        
        % Flip the axes if onset response is shown
        flag_flip_onset_axes = 0;
        if (flag_flip_onset_axes && contains(figure_title_str,'onset'))
            set(gca,'YDir','reverse')
            set(gca,'XDir','reverse')
            set(gca,'XAxisLocation','top')
            set(gca,'YAxisLocation','right')
        end
        
        %     figure; histogram(x-y)
        
    end
    
    