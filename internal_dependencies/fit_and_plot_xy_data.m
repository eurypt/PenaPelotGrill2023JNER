%{


Guide to inputs:
flag_scale_std_by_y:
set to false (default) or true;  if set to false, the std of the fit will
be plotted as constant across the entire domain; if set to true, the std of
the fit will be plotted as dependent on the y value

std_type:
string that can be set to either 'sample' (default) or 'mean'; if set to
'sample', the std of individual new samples from the fit will be plotted;
if set to mean, the std of the mean of many samples from the fit will be
plotted

flag_plot_scatter_data:
set to false or true (default);  if set to true, the raw data points (x &
y) will be plotted as scatter against the fit

fig_handle:
variable containing a figure handle; if specified, the plot will be
performed on that figure handle rather than creating a new figure window

line_color:
RGB color vector with values from 0 to 1 (e.g., [0 0 0]) indicating the
color to make the line; if not specified, by default the colors will be
based on MATLAB 2018a default line colors
%}

function fit_obj = fit_and_plot_xy_data(x,y,x_shift,fit_obj,gof,std_type,flag_scale_std_by_y,flag_plot_scatter_data,fig_handle,line_color)

%%% Check inputs
if (~exist('std_type','var'))
    std_type = 'sample';
end
assert(strcmp(std_type,'sample') || strcmp(std_type,'mean'),'invalid value of std_type')

if (~exist('flag_scale_std_by_y','var'))
    flag_scale_std_by_y = false;
end
assert(flag_scale_std_by_y==0 || flag_scale_std_by_y==1 || flag_scale_std_by_y==true || flag_scale_std_by_y==false,...
    'Invalid value of flag_scale_std_by_y');

if (~exist('flag_plot_scatter_data','var'))
    flag_plot_scatter_data = true;
end
assert(flag_plot_scatter_data==0 || flag_plot_scatter_data==1 || flag_plot_scatter_data==true || flag_plot_scatter_data==false,...
    'Invalid value of flag_plot_scatter_data');

if (~exist('fig_handle','var') || isempty(fig_handle))
    fig_handle = figure('color',[1 1 1],'position',[680,517,631,461]);
else
    figure(fig_handle); % set fig_handle as the current active figure
end

if (~exist('line_color','var') || isempty(line_color))
    % determine number of lines already present in the figure
    gg=get(gca,'Children');
    number_of_lines_present = length(gg);
    
    % pick color based on MATLAB 2018a default line colors
    MATLAB_2018a_default_colors = [
        0    0.4470    0.7410
        0.8500    0.3250    0.0980
        0.9290    0.6940    0.1250
        0.4940    0.1840    0.5560
        0.4660    0.6740    0.1880
        0.3010    0.7450    0.9330
        0.6350    0.0780    0.1840
        ];
    line_color = MATLAB_2018a_default_colors(mod(number_of_lines_present,size(MATLAB_2018a_default_colors,1))+1,:);
end

fit_coeffs = coeffvalues(fit_obj)';

std_y = gof.rmse; % Same as the following manual process: SSE = sum((ones(length(y),1)-[((x-x_shift).^2)./y, (x-x_shift)./y, 1./y]*fit_coeffs).^2); MS_Re = SSE/(length(y)-3); std_y = sqrt(MS_Re);

%%% Make x & y values & std vs. x values for output so user can interpolate
%%% the appropriate y and std values for a given x value
x_fit_vals = (linspace(min(x),max(x),1000))'; % make interpolant values from 1.12 um up to 20 um, since below this value the geometry is invalid anyway due to negative STIN lengths; above that, there aren't any fibers anyway
y_fit_vals = [(x_fit_vals-x_shift).^2, (x_fit_vals-x_shift), ones(length(x_fit_vals),1)]*fit_coeffs;
SSX = sum((x-mean(x)).^2);


if (strcmp(std_type,'sample'))
    std_mean_fit_eval = std_y*sqrt(1 + 1/length(x) + ((x_fit_vals-x_shift).^2)./SSX);
elseif (strcmp(std_type,'mean'))
    std_mean_fit_eval = std_y*sqrt(1/length(x) + ((x_fit_vals-x_shift).^2)./SSX);
end

if (flag_scale_std_by_y)
    std_mean_fit_eval = y_fit_vals.*std_mean_fit_eval;   
end


% subaxis(1,1,1,'MT',0.17,'MB',0.2,'ML',0.2)

% (variance of the mean of fit) = (variance about mean) + (variance due to distance from mean x)
df = length(x)-3;
t_value_span = linspace(0,5,1000); p_value_span = 1-abs(tcdf(t_value_span,df)-tcdf(-t_value_span,df));
alpha = 0.05;
[~,p_idx] = min(abs(p_value_span-alpha));
ci_eval_positive = std_mean_fit_eval*t_value_span(p_idx);
ci_eval_negative = ci_eval_positive;
shadedErrorBar(x_fit_vals,y_fit_vals,[ci_eval_positive';ci_eval_negative'],...
    'lineprops',{'-','color',line_color},'transparent',true) % plot using confidence intervals of beta
hold on
if (flag_plot_scatter_data)
    plot(x,y,'k.','MarkerSize',15)
end
x_offset = 1; % [um]
xlim([min(x_fit_vals)-x_offset,max(x_fit_vals)+x_offset])
set(gca,'FontSize',14)



%{
Incomplete / Archived / Obsolete 

%{


string that can be set to either 'flat' (default) or 'yscaled'; if set to
'flat', the std of the fit will be plotted as constant across the entire
domain; if set to 'yscaled', the std of the fit will be plotted as
dependent on the y value


%%% Check inputs
if (~exist('flag_scale_std_by_y','var'))
    flag_scale_std_by_y = 'flat';
end
assert(strcmp(flag_scale_std_by_y,'flat') || strcmp(flag_scale_std_by_y,'yscaled'),'invalid value of flag_scale_std_by_y')

if (strcmp(flag_scale_std_by_y,'yscaled'))
    std_mean_fit_eval = y_fit_vals.*std_mean_fit_eval;   
end

%}

%{
std_scaling:
a vector of the same size as y that gets scaled by the std of the fit; by
default, std_scaling is just a vector of ones; specifying this vector
enables the code to plot error bars in the case where uncertainty depends
on the value of the measurements
%}

%{
%%% Check inputs
if (~exist('std_scaling','var'))
    std_scaling = ones(size(y));
end
assert(isequal(size(std_scaling,y)))
%}

%{
%%% Plot equation itself
fit_equation_str1 = ['$y = (',...
    num2str(fit_coeffs(1),'%0.2f') '\frac{1}{\mu m} (x-',num2str(x_shift,'%0.2f'),' \mu m)^2 + $'];
fit_equation_str2 = ['$ ',num2str(fit_coeffs(2),'%0.2f'), ' (x-',num2str(x_shift,'%0.2f'),' \mu m) + ',...
    num2str(fit_coeffs(3),'%0.2f'), ' \mu m)*N(1,',num2str(std_y,'%0.2f'),')$'];

font_size = 16;
text(-0.086, 1.1913,fit_equation_str1,'units','normalized','FontSize',font_size,'interpreter','latex')
text(-0.0057, 1.0789,fit_equation_str2,'units','normalized','FontSize',font_size,'interpreter','latex')
set(gca,'FontSize',16)
%}

%}