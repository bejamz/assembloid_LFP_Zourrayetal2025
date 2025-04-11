function [inst_amp,inst_phase] = processHilberts(data,si_us)

%%% data should be bandpass filtered to the frequency of interest 
%%% (for biological data; not if single dominant frequency exists throughout)
%%% see https://www.gaussianwaves.com/2017/04/extract-envelope-instantaneous-phase-frequency-hilbert-transform/
%%%
%%% Inputs:
%%%     - data - nx1 vector
%%%     - si_us - sampling interal IN MICROSECONDS!
%%%
%%% FUN FACT: hilbert(x) does not compute the hilbert transform of x...

fprintf('... running Hilbert analyses\n')

z = hilbert(data); % Extract analytical signal
inst_amp = abs(z);
inst_phase = mod(unwrap(angle(z)),2*pi);
inst_freq = diff(inst_phase)/(2*pi)*(1000000/si_us);%inst frequency
