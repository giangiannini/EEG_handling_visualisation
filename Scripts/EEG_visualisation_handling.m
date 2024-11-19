%% Add the fieldtrip path
addpath('C:/Users/nnu02/Documents/MATLAB/fieldtrip-20220104'); %depends on your Fieldtrip Path.

folder = 'C:\Users\<yourUsernamehere>\Desktop\EEG_handling_visualisation\'; %Windows
folder = '/Users/<yourUsernamehere>/Desktop';
dataset = strcat(folder, 'ID01\ID01.bdf');

%% Load the raw dataset
cfg = []; 
cfg.continuous = 'yes'; 
cfg.dataset = dataset;
continuous = ft_preprocessing(cfg); 

% Downsample so it's more manageable
cfg = []; 
cfg.resamplefs = 256; 
continuous = ft_resampledata(cfg, continuous); 

%% Open the raw dataset
cfg = [];
cfg.continuous = 'yes';
cfg.ylim = [-20 20];
cfg.blocksize = 5;
%cfg.datafile = strcat(folder, '24_exp_SPN-Cropped.bdf');
ft_databrowser(cfg, continuous);

%% Try to adjust it a bit :) 
cfg = []; 
cfg.reref = 'yes'; 
cfg.refchannel = 'EEG'; 
cfg.hpfilter = 'yes'; %HPF
cfg.hpfreq = 0.1; %HPF = 0.1 Hz
cfg.hpfilttype = 'but'; %Butterworth
cfg.hpfiltdir = 'twopass'; %Two-pass filter (to avoid distortions in both directions)
cfg.hpinstabilityfix = 'reduce'; %Filter order is determined automatically based on instability.
cfg.padding = size(continuous.trial{1,1},2)/continuous.fsample +100*2; %We also need extra padding
continuous_reref_hpf = ft_preprocessing(cfg, continuous); 

cfg = [];
cfg.continuous = 'yes';
cfg.ylim = [-20 20];
cfg.blocksize = 5;
cfg.channel = 1:64;
ft_databrowser(cfg, continuous_reref_hpf);
%% Epoching!
%First let's have a look at our events of interest
cfg                         = [];
cfg.dataset                 = dataset;
cfg.trialfun                = 'ft_trialfun_general'; % this is the default
cfg.trialdef.eventtype      = 'STATUS';
cfg.trialdef.eventvalue     = 128; % the values of the stimulus trigger for the three conditions
prova = ft_read_event(dataset);
for i = 1:length(prova)
    prova(i).sample = round(prova(i).sample/8);
end

continuous_reref_hpf.cfg.previous = rmfield(continuous_reref_hpf.cfg.previous, 'origfs');

cfg = []; 
cfg.event = prova;
cfg.plotevents = 'yes'; 
cfg.continuous = 'yes';
cfg.ylim = [-20 20];
cfg.blocksize = 5;
cfg.channel = 1:64;
ft_databrowser(cfg, continuous_reref_hpf);

%Then do proper epoch
trl = [];
cfg                         = [];
cfg.dataset                 = dataset;
cfg.trialfun                = 'ft_trialfun_general'; % this is the default
cfg.trialdef.eventtype      = 'STATUS';
cfg.trialdef.eventvalue     = 128; % the values of the stimulus trigger for the three conditions
cfg.trialdef.prestim        = 0.1; % in seconds
cfg.trialdef.poststim       = 0.5; % in seconds
cfg = ft_definetrial(cfg);
trl = cfg.trl; 

cfg = [];
cfg.trl = trl;
cfg.trl(:,1:3) = round(cfg.trl(:,1:3)/8);
epoched = ft_redefinetrial(cfg, continuous_reref_hpf);

cfg = []; 
cfg.ylim = [-20 20];
cfg.channel = 1:64;
ft_databrowser(cfg, epoched);

%% TRIAL REJECTION
load("visual.mat")
cfg = []; 
cfg.ylim = [-20 20];
cfg.channel = 1:64;
cfg.artfctdef.visual.artifact = visual; 
artifact = ft_databrowser(cfg, epoched);
visual = artifact.artfctdef.visual.artifact;
save("visual.mat", "visual"); %save and load every time we run this section.

cfg = []; 
cfg.artfctdef.reject = 'complete'; %complete removal of trials contaminated
cfg.artfctdef.visual.artifact = visual;
cleaned_epoched = ft_rejectartifact(cfg, epoched); 

%% AVERAGE OVER TRIALS
cfg = []; 
cfg.channel = 1:64; 
ERP = ft_timelockanalysis(cfg, cleaned_epoched);

cfg = []; 
cfg.ylim = [-0.5 0.5];
ft_databrowser(cfg, ERP);

cfg = []; 
cfg.ylim = [-8 8];
cfg.viewmode = 'butterfly';
ft_databrowser(cfg, ERP);

%% BASELINE 
cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow = [-0.1 0];
ERP_baseline = ft_preprocessing(cfg, ERP);

cfg = [];
cfg.ylim = [-8 8];
cfg.viewmode = 'butterfly';
ft_databrowser(cfg, ERP_baseline);

%% EXPLORATION OF THE RESULTS 
caplocation = strcat(folder, 'biosemi64.lay');

cfg = [];
cfg.layout = caplocation;
ft_multiplotER(cfg, ERP_baseline)




















%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% FOR PLOTTING OTHER IMAGES IN THE PRESENTATION %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATE DATA (simple)
time = -1:0.0001:1; 

A = 2; 
sin10 = A*sin(2*pi*10*time); 
sin5 = A*sin(2*pi*5*time); 
sin30 = A*sin(2*pi*30*time); 
total_EEG = sin10 + sin5 + sin30; 

figure;
subplot(4,1,1)
plot(time, sin5)
title('Sine wave of amplitude 2 and 5Hz')
subplot(4,1,2)
plot(time, sin10)
title('Sine wave of amplitude 2 and 10Hz')
subplot(4,1,3)
plot(time, sin30)
title('Sine wave of amplitude 2 and 30Hz')
subplot(4,1,4)
plot(time, total_EEG); 
title('Linear summation of the signals')
xlabel('Time')
ylabel('Voltage')

N = numel(time);
for i = 1:40
    sine_wave = A*sin(2*pi*(i-1)*time); 
    fourier(i) = sum(sine_wave.*total_EEG); 
end
fourier = fourier / N;
plot((1:40)-1, fourier)
xlabel('Frequency')
ylabel('Power')
title('Fourier transform on the "Total EEG"')

%% GENERATE MORE REALISTIC

%Resampling / Nyquist
time = -1:0.0001:1; 
A = 2;
figure; 
plot(time, A*sin(2*pi*2*time)); 
hold on
time = linspace(-1,1,201); %sampling 100Hz
plot(time, A*sin(2*pi*2*time), 'r*')
time = linspace(-1, 1, 5); %sampling 4Hz
plot(time, A*sin(2*pi*2*time), 'b*')
legend({'Wave at 2Hz', 'Sampling at 100Hz', 'Sampling at 2Hz'})


%1/f plot
frequency = 1:0.01:50;
f = 1./frequency; 
figure; 
plot(frequency, f); 
ylabel('power')
xlabel('frequency')
title('Pink noise example')

A = 2; 
sin10 = A/10*sin(2*pi*10*time); 
sin5 = A/5*sin(2*pi*5*time); 
sin30 = A/30*sin(2*pi*30*time); 
total_EEG = sin10 + sin5 + sin30; 

figure;
subplot(4,1,1)
plot(time, sin5)
title('Sine wave of amplitude 2/5 and 5Hz')
subplot(4,1,2)
plot(time, sin10)
title('Sine wave of amplitude 2/10 and 10Hz')
subplot(4,1,3)
plot(time, sin30)
title('Sine wave of amplitude 2/30 and 30Hz')
subplot(4,1,4)
plot(time, total_EEG); 
title('Linear summation of the signals')
xlabel('Time')
ylabel('Voltage')


figure;
N = numel(time);
for i = 1:40
    sine_wave = A*sin(2*pi*(i-1)*time); 
    fourier(i) = sum(sine_wave.*total_EEG); 
end
fourier = fourier / N;
plot((1:40)-1, fourier)
xlabel('Frequency')
ylabel('Power')
title('Fourier transform on the "Total EEG", adjusted for 1/f')

%% GENERATE EVEN MORE REALISTIC
A = 2; 
sin01 = A/0.001*(-1)*sin(2*pi*0.001*time);
sin10 = A/10*sin(2*pi*10*time); 
sin5 = A/5*sin(2*pi*5*time); 
sin30 = A/30*sin(2*pi*30*time); 
total_EEG = sin01 + sin10 + sin5 + sin30; 

figure;
subplot(5,1,1)
plot(time, sin01); 
title('Sine wave of amplitude 2/0.01 and 0.01 Hz')
subplot(5,1,2)
plot(time, sin5)
title('Sine wave of amplitude 2/5 and 5Hz')
subplot(5,1,3)
plot(time, sin10)
title('Sine wave of amplitude 2/10 and 10Hz')
subplot(5,1,4)
plot(time, sin30)
title('Sine wave of amplitude 2/30 and 30Hz')
subplot(5,1,5)
plot(time, total_EEG); 
title('Linear summation of the signals')
xlabel('Time')
ylabel('Voltage')

fourier = []; 
N = numel(time);
frequencies = 0:0.1:40;
for i = 1:length(frequencies)
    sine_wave = A*sin(2*pi*(frequencies(i))*time); 
    fourier(i) = sum(sine_wave.*total_EEG); 
end
figure;
fourier = fourier / N;
plot(frequencies(1:10:length(frequencies)),abs(fourier(1:10:length(fourier))))
xlabel('Frequency')
ylabel('Power')
title('Fourier transform on the "Total EEG", plus slow freqs')

time = -1:0.01:1; %Let's downsample the signal to 100Hz so it's manageable for simulations
A = 3; 
signal = A*sin(2*pi*1*time); 
noise = rand(1,1)*A*sin(2*pi*rand(1,1)*10*time); 
total_noise = []; 
for i = 1:100 %change to 10 and 100 and see the difference.
    noise = A*sin(2*pi*rand(1,1)*5*time + 100*rand(1,1)); 
    noise = signal + noise; 
    total_noise = [total_noise; noise]; 
end

averaged_signal_noise = mean(total_noise, 1);
figure; 
subplot(2,1,1)
plot(time, signal);
title('original signal')
subplot(2,1,2)
plot(time, averaged_signal_noise); 
title('original signal + random noise (100 trials)')
xlabel('time')
ylabel('voltage')

averaged_signal_noise = [averaged_signal_noise(1:100)-signal(1:100) averaged_signal_noise(101:201)*1];

figure; 
plot(time, averaged_signal_noise); hold on
plot(time, averaged_signal_noise+3); 
legend({'Baseline corrected', 'no Baseline correction'})









