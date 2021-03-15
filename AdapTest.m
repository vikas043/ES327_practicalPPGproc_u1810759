function PPGsig= AdapTest()
clear all;
%% Part 1-Generating a known signal %%
fs=256;
t=0:1/fs:40;
f=6;
FFTres=1024;
WFlength = 16; % Wiener filter length

sinePPG=4*sin(2*pi*f*t);

%% Part 2- Simulating accelerometer noise %%

Acc_x=(2.2*randn(size(sinePPG)))+0.05;
Acc_y=(0.9*randn(size(sinePPG)));
Acc_z=(1.8*randn(size(sinePPG)))+0.02;

%% Part 3- Viewing a noisy signal %%
%Not actually representative for MA,not used for filterimplementation
%Just a visualizer
sinePPGnoisy=sinePPG+1/3*(Acc_x+Acc_y+Acc_z);
figure()
plot(t,sinePPG);grid on;
hold on;
plot(t,sinePPGnoisy);
title('Pure and noisy input signal');
xlabel('Time (s)');
ylabel('PPG signal [-]');
xlim([0 3]);
legend({'Pure input sinusoid','Noisy input sinusoid'})
figure()
pspectrum(sinePPG);
figure()
pspectrum(sinePPGnoisy);

sinePPGnoisyStd=(sinePPGnoisy-mean(sinePPGnoisy))/(std(sinePPGnoisy));
% framework rule, 4s window 2s shift
window   = 4 * fs;  % window length is 4 seconds
step     = 1 * fs;  % step size is 2 seconds


wLength = floor((length(sinePPGnoisy)-window)/step) + 1;  % total number of windows(estimates)
    
    clear W11_FFTi W11_PPG_ave_FFT_Clean;
    
    j=1;
    for i =  [1 :  wLength]
        curSeg = (i-1)*step+1 : (i-1)*step+window;
        PPG_ave = sinePPGnoisyStd(:,curSeg); PPG_ave2 = sinePPGnoisy(:,curSeg);
        ACC_X_seg = Acc_x(:,curSeg); ACC_Y_seg = Acc_y(:,curSeg); ACC_Z_seg = Acc_z(:,curSeg);
        

        % PPG FFT
        PPG_ave_FFT = fft(PPG_ave,FFTres);
        PPG_ave2_FFT = fft(PPG_ave2,FFTres);
        
        %acceleration FFT
        ACC_X_FFT= fft(ACC_X_seg,FFTres);
        ACC_Y_FFT= fft(ACC_Y_seg,FFTres);
        ACC_Z_FFT= fft(ACC_Z_seg,FFTres);
        
        %Reference A Temko "Accurate Heart Rate Monitoring During Physical Exercises Using PPG," 
        %in IEEE Transactions on Biomedical Engineering,vol. 64,no. 9,pp. 2016-2024,Sept. 2017      
        % Wiener filtering PPG-ACC, two types
        % W1 filter = clean_signal/(clean_signal+noise) = (all_signal-noise)/all_signal,
        % where all_signal = time_average(ppg_signal), noise = time_average(accelerometer)
        %
        % The process of Wiener filtering:
        % 1. Estimate the 'signal'/'clean_signal' level as average over
        %    past 15 spectral envelopes (~30s)
        % 2. Scale the spectrum of noise and signal/clean_signal
        % 3. Apply the filter as weighting of the spectral components
          
        
%Normalized filtering for HR estimation
        % Algorithm 1, spectral envelope calculated using second power of FFT
        %adding power spectrum of current window signal to accumalate input spectrum
        W11_FFTi(i,:) = abs(PPG_ave_FFT).^2;
        
        %averaged input power spectrum
        %edge case for the first window
        if i==1, W11_PPG_ave_FFT_ALL = W11_FFTi(i,:); 
        %averaged input power spectrum using upto 16 windows
        else W11_PPG_ave_FFT_ALL = mean(W11_FFTi(max(1,i-WFlength):i,:),1); 
        end
        
        %normalizing power spectrums
        W11_PPG_ave_FFT_ALL_norm = (W11_PPG_ave_FFT_ALL)/max(W11_PPG_ave_FFT_ALL);
        W11_ACC_X_FFT_norm = (abs(ACC_X_FFT).^2)/max(abs(ACC_X_FFT.^2));
        W11_ACC_Y_FFT_norm = (abs(ACC_Y_FFT).^2)/max(abs(ACC_Y_FFT).^2);
        W11_ACC_Z_FFT_norm = (abs(ACC_Z_FFT).^2)/max(abs(ACC_Z_FFT).^2);
        %Average power spectrum of accelerometer signals
        W11_ACC_N_FFT_norm=1/3*(W11_ACC_X_FFT_norm+W11_ACC_Y_FFT_norm+W11_ACC_Z_FFT_norm);
        %Filter weights spectrum for the particular window 
        WF11 = (1 - (W11_ACC_N_FFT_norm./(W11_PPG_ave_FFT_ALL_norm)));
        %Adding Clean signal estimate to the output
        PPG_ave_FFT_Clean(i,:) = abs(PPG_ave_FFT).*(WF11);

        
 %Unnormalized filtering for output plotting and SpO2
 
 %adding power spectrum of current window signal to accumalate input spectrum
        W12_FFTi(i,:) = abs(PPG_ave2_FFT).^2;
 %edge case for the first window       
        if i==1, W12_PPG_ave_FFT_ALL = W12_FFTi(i,:);
 %averaged input power spectrum using upto 16 windows
        else W12_PPG_ave_FFT_ALL = mean(W12_FFTi(max(1,i-WFlength):i,:),1); 
        end
        %Average power spectrum of accelerometer signals
        W12_ACC_N_FFT=1/3*((abs(ACC_X_FFT).^2)+(abs(ACC_Y_FFT).^2)+(abs(ACC_Z_FFT).^2));
        %Filter weights spectrum for the particular window 
        WF12 = (1 - (W12_ACC_N_FFT./(W12_PPG_ave_FFT_ALL)));
        %Adding Clean signal estimate to the output
        PPG_ave2_FFT_Clean(i,:) = abs(PPG_ave2_FFT).*(WF12);
        
        
        if mod(i-1,4)==0
        PPG_ave_IFFT_Clean(j,:)= ifft(PPG_ave_FFT.*(WF11));
        PPG_ave2_IFFT_Clean(j,:) = ifft(PPG_ave2_FFT.*(WF12));
        j=j+1;
        end
    end
    
%reshaping normalized clean PPG to a vector
PPG_Clean_normt=PPG_ave_IFFT_Clean.';
PPG_Clean_norm=reshape(PPG_Clean_normt,1,[]);
%reshaping clean PPG to a vector
PPG_Cleant=PPG_ave2_IFFT_Clean.';
PPG_Clean=reshape(PPG_Cleant,1,[]);

t2=0:1/fs:(length(PPG_Clean)/fs);%for plotting, time values from sample rate and index
t2=t2(2:length(t2));
%Plotting figures of clean PPG
figure();
plot(t2,PPG_Clean_norm);grid on;
title('Processed signal after normalization');
xlabel('Time (s)');
ylabel('PPG signal [-]');
xlim([0 3]);
figure();
plot(t2,PPG_Clean);grid on;
title('Processed signal');
xlabel('Time (s)');
ylabel('PPG signal [-]');
xlim([0 3]);

PPGsig=PPG_Clean;

