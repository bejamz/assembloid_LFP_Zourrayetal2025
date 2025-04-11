function plotOrgSummaryFigure(organoid,savefile)

h = figure;
ch = 1;
baseline_min = 15;

data = organoid.Data.data(:,1);

%% First let's plot the trace itself

ax = subplot(4,3,[1:3]);


time = [1:length(data)]'./(organoid.Data.Fs*60);
pl = plot(ax,time,data);

xlim([0 120])
ylim([nanmean(data)-20*nanstd(data), nanmean(data)+20*nanstd(data)])
ylabel('Voltage /mV');
xlabel('Time /min'); 

if isfield(organoid.Meta,'Group')
if organoid.Meta.Group == 1
    pl.Color = 'r';
elseif organoid.Meta.Group == 2
    pl.Color = 'k';
end
end

%% Let's plot the PSD

baseline_ind = round(baseline_min*60*organoid.Data.Fs);

[pxx,fxx] = pwelch(data(1:baseline_ind),64000,32000,[],organoid.Data.Fs); % somewhat arbitrary parameters
[pxx2,fxx2] = pwelch(data(baseline_ind+1:end),64000,32000,[],organoid.Data.Fs); % somewhat arbitrary parameters

ax = subplot(4,3,4);
hold on;
title('Power spectral densities')
plot(fxx,pxx);
plot(fxx2,pxx2);
xlim([0 50]);
legend('Baseline PSD','+55mM KCl PSD')

ax = subplot(4,3,7);
hold on;
plot(fxx,pxx);
plot(fxx2,pxx2);
xlim([0 50]);
ax.YScale = 'log';


ax = subplot(4,3,10);
freqs = [1:2:101];
psd_fold_change = NaN(length(freqs)-1,1);
for ind = 1:length(freqs)-1
    psd_fold_change(ind) = nanmean(pxx2(fxx2>=freqs(ind)&fxx2<freqs(ind+1))) ./ nanmean(pxx(fxx>=freqs(ind)&fxx<freqs(ind+1)));
end
plot(ax,freqs(1:end-1),psd_fold_change);

%%

ax = subplot(4,3,6);
if isfield(organoid.SmallBlockAnalysis,'small_block_analysis_done') && organoid.SmallBlockAnalysis.small_block_analysis_done
    plot(organoid.SmallBlockAnalysis.stdevs_norm);
    title('STDev in small blocks')
else
    text(0.3,0.5,'No perturbation/normalised analysis')
    axis off
end

ax = subplot(4,3,9);
plot(organoid.SmallBlockAnalysis.norm_coastline_smallblocks);
title('Norm coastline in small blocks')

ax = subplot(4,3,5);
hold on;
plot(organoid.LargeBlockAnalysis.bandpow_theta);
plot(organoid.LargeBlockAnalysis.bandpow_sgamma);
plot(organoid.LargeBlockAnalysis.bandpow_fgamma);
plot(organoid.LargeBlockAnalysis.bandpow_hfos);
legend('Theta power','S Gamma power','F gamma power','HFO power')
title('Raw power in large blocks')

ax = subplot(4,3,8);
hold on;
plot(organoid.LargeBlockAnalysis.proportion_theta);
plot(organoid.LargeBlockAnalysis.proportion_sgamma);
plot(organoid.LargeBlockAnalysis.proportion_fgamma);
plot(organoid.LargeBlockAnalysis.proportion_hfos);
title('Proportional power in large blocks')

if exist('savefile','var')
    
    h.Units = 'normalized';
    h.Position = [0,0,1,1];
    saveas(h,savefile);

end

