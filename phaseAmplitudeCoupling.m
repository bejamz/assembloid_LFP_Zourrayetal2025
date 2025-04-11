function [phaseamps,phase_bins] = phaseAmplitudeCoupling(phases,amps)

assert(length(phases)==length(amps))

phase_bins = [0:1:359]+0.5;

phases = round(rad2deg(phases));
phases = mod(phases,360);

phaseamps = NaN(360,1);
for th = 0:359;
    phaseamps(th+1) = nanmean(amps(phases==th));
end

plot(phaseamps)