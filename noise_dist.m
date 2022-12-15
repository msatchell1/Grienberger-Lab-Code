function [S_all] = noise_dist(S_all, view_fits)
% Seperates the base noise of dF/F signals from activity by plotting dF/F
% points on a histogram, than calculates the standard deviation of this
% noise signal. Does so independently for each neuron.
%
% Also takes input view_fits. 1 to toggle each gaussian fit, 0 to skip
% through to the last fit.
%
% ----- Michael Satchell 12/14/22 -------

% Use the smoothed dF/F data:
dFF_data = S_all.datasetSm;

num_nrns = size(dFF_data,2);

% Creates array of zeros to hold the standard deviation of noise signal for 
% each neuron.
S_all.dFF_noise_std = zeros(1,num_nrns);

fig1 = figure(1);

for i = 1:num_nrns
    
    nrn_dFF = dFF_data(:,i); % Fluorescence data for a single neuron.
    
    clf(fig1);
    hold on;
    H = histogram(nrn_dFF, 500); % Plotting the data with 500 bins.

    xvals_hist = H.BinEdges(1:end-1)'; % Approximating the left edges of each 
    % bin as the location of the bin.
    yvals_hist = H.Values'; % Number of counts in each bin.

    [max_val, max_i] = max(yvals_hist); % Finds the largest bin to center the gaussian around.

    fo = fitoptions('Method','NonlinearLeastSquares',... % Squared error fit.
      'Lower',[-1,0,-inf,-inf],...   % Lower bounds for [a b c d].
      'Upper',[1,inf,inf,inf],... % Upper bounds for [a b c d].
      'Startpoint',[0, max_val, xvals_hist(max_i), 0.1]); % Startpoint for the fit. 
    % There should be no vertical offset a, the amplitude b should be close to 
    % the maximum bin value, the centroid of the curve c should be near the 
    % location of the largest bin, and the standard deviation d should be somewhere
    % on the order of 0.1.
    
    f1 = fittype('a+b*exp(-((x-c).^2)/((2*d^2)))', 'options',fo); % Loading 
    % the fit options into the
    
    [gausscurve, gofgauss] = fit(xvals_hist, yvals_hist, f1); % Performs fit

    % Adds standard deviation calculated from fit to array.
    S_all.dFF_noise_std(1,i) = gausscurve.d;
    
    % Plots each fit individually with R-squared value.
    yfitdata = gausscurve(min(xvals_hist):H.BinWidth:max(xvals_hist));
    p = plot(min(xvals_hist):H.BinWidth:max(xvals_hist), yfitdata);
    p.LineWidth = 1;
    ylabel("Bin Counts")
    xlabel("dF/F")
    title(strcat('dF/F Histogram and Gauss Fit for Cell', num2str(i)));
    % Add r-squared value to plot.
    txt = strcat('R^2 = ', num2str(round(gofgauss.rsquare,3)));
    text(max(xvals_hist)/2,max(yvals_hist),txt);
    
    if view_fits
        pause; % Use pause to view each plot.
    end

end



end