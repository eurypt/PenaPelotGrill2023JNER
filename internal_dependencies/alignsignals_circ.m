%{

Align periodic signals based on maximizing cross correlation.
Similar functionality to the built-in function alignsignals, except that
this is specifically designed for periodic signals

top_N_waveforms - a N-by-T matrix, where N is the number of signals and T
is the time component.

%}

function top_N_waveforms_aligned = alignsignals_circ(top_N_waveforms)

align_method = 2;
top_N_waveforms_aligned = zeros(size(top_N_waveforms));
for j = 1:size(top_N_waveforms_aligned,1)
    switch align_method
        case 1
            [~,top_N_waveforms_aligned(j,:)] = alignsignals(top_N_waveforms(1,:),top_N_waveforms(j,:));
        case 2
            % Calculate all circularly shifted version of the
            % non-reference waveform
            n_timepoints = length(top_N_waveforms(j,:));
            lags = 0:(n_timepoints-1);
            non_reference_waveform_circshifted_matrix = cell2mat(arrayfun(@(x) ...
                circshift(top_N_waveforms(j,:),x)', lags,'UniformOutput',false))';
            
            % Calculate xcorr manually between reference waveforms and
            % non-reference circularly shifted waveforms
            flag_make_reference_waveform_sine = 0;
            if (flag_make_reference_waveform_sine)
                reference_waveform = cos(linspace(0,2*pi-(2*pi/size(top_N_waveforms,2)),size(top_N_waveforms,2))-pi/2);
            else
                reference_waveform = top_N_waveforms(1,:);
            end
            xcorr_values = reference_waveform*(non_reference_waveform_circshifted_matrix');
            
            % Store the circularly shifted non-reference waveform with
            % the largest xcorr value
            [~,max_idx] = max(xcorr_values);
            top_N_waveforms_aligned(j,:) = non_reference_waveform_circshifted_matrix(max_idx,:);
            
        otherwise
            error('Invalid align_method')
    end
end