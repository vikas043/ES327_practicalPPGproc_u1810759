 function PPGsigOut= WienerPPGfirMar()
clear all;
%% Part 1-Setting up raw signals %%
%%Comment to change the subject 

%%Load files containing PPG sensor readings,ECG and accelerometer data
raw=load('S3_step1_ppg.mat');
raw_acc=load('S3_step1_acc.mat');

% Uncomment to change the activity
% 
% raw=load('S3_rest1_ppg.mat');
% raw_acc=load('S3_rest1_acc.mat');
%
%Uncomment to change the subject
%
% raw=load('S6_squat2_ppg.mat');
% raw_acc=load('S6_squat2_acc.mat');
% 
% Uncomment to change the activity
% 
% raw=load('S6_rest2_ppg.mat');
% raw_acc=load('S6_rest2_acc.mat')
% 

% Comment when using resting data
% Extracting data in to vectors
rawPPG=raw.PPG(:,2);%PPG signal from PPG file
rawPPG=rawPPG.';
%ACC signals from ACC file
ACC_Ref_c3=raw_acc.ACC(:,2);
ACC_Ref_c3=ACC_Ref_c3.';
ACC_Ref_c4=raw_acc.ACC(:,3);
ACC_Ref_c4=ACC_Ref_c4.';
ACC_Ref_c5=raw_acc.ACC(:,4);
ACC_Ref_c5=ACC_Ref_c5.';

% Uncomment when using resting data
% rawPPG=raw.PPG(1:24000,2);%PPG signal from PPG file
% rawPPG=rawPPG.';
% %ACC signals from ACC file
% ACC_Ref_c3=raw_acc.ACC(1:24000,2);
% ACC_Ref_c3=ACC_Ref_c3.';
% ACC_Ref_c4=raw_acc.ACC(1:24000,3);
% ACC_Ref_c4=ACC_Ref_c4.';
% ACC_Ref_c5=raw_acc.ACC(1:24000,4);
% ACC_Ref_c5=ACC_Ref_c5.';

%Standardized PPG
rawPPGav1=(rawPPG-mean(rawPPG))/(std(rawPPG));
%Untandardized PPG
rawPPGav2=rawPPG;

fs=400;%sampling frequency
t=0:1/fs:(length(rawPPG)/fs);%for plotting, time values from sample rate and index
t=t(2:length(t));

%Plots after first stage of preprocessing- averaging different channels
figure();
plot(t,rawPPGav1);grid on;
title({'Standardised Raw PPG signal','Subject 1'});
xlabel('Time (s)');
ylabel('PPG signal[-]');
figure();
plot(t,rawPPGav2);grid on;
title({'Raw PPG signal','Subject 1'});
xlabel('Time (s)');
ylabel('PPG signal[na]');


%% Part 2-Detrending the signals %%
detPPGav1=detrend(rawPPGav1);
detPPGav2=detrend(rawPPGav2);
figure();
plot(t,detPPGav1);grid on;
title({'Standardized Detrended PPG signal','Subject 1'});
xlabel('Time (s)');
ylabel('PPG signal[-]');
figure();
plot(t,detPPGav2);grid on;
title({'Detrended PPG signal','Subject 1'});
xlabel('Time (s)');
ylabel('PPG signal[-]');


%% Part 3-Butterworth Filter %%
fn=fs/2;%Nyquist frequency
fclow=0.5/fn;%PPG lower cut-off frequency
fchigh=4.5/fn;%PPG higher cut-off frequency

%Array of cut-off frequencies
Wn=[fclow fchigh];

%Filter order trial (same as butterworth)
n=8;
%generating filter coefficients for Hamming window filter
b=fir1(n,Wn);
%Plotting filter frequency response
figure();
freqz(b,1);
%Recommended filter order for FIR bandpass filter
n=64;
%generating filter coefficients
b=fir1(n,Wn);
%Plotting filter frequency response
figure();
freqz(b,1);
%filtering the PPG signals
filPPGav1=filter(b,1,detPPGav1);
filPPGav2=filter(b,1,detPPGav2);

%Preprocessing accelerometer data for Adaptive stage
ACC_X = filter(b,1,ACC_Ref_c3);
ACC_Y = filter(b,1,ACC_Ref_c4);
ACC_Z = filter(b,1,ACC_Ref_c5);
%writing to output
PPGsigOut.filPPGstdz=filPPGav1;
PPGsigOut.filPPG=filPPGav2;

%% Part 4-BPM calculations %%
%saves the amplitude and index of peaks in separate vectors
%filters peaks and reduces False Positives to the extent of sanity in
%expected parameters, does not unfairly reduce HR, this can be verified
%visually and on the basis of sampling rate

[peaks,locs]=findpeaks(filPPGav1,'MinPeakHeight',0.1,'MinPeakDistance',90,'MinPeakProminence',0.03);
%time at peak locations
locs_t=locs/fs;
%Plotting after second stage of preprocessing and bandpass filtering
%Plotting detected peaks
figure();
plot(t,filPPGav1);
hold on;
plot(locs_t,filPPGav1(locs),'ro','MarkerFaceColor','red','MarkerSize',3);
grid on;
title({'Standardized filtered PPG signal','Subject 1'});
xlabel('Time (s)'); 
ylabel('PPG signal[-]');


%saves amplitude and index of peaks
[peaks2,locs2]=findpeaks(filPPGav2,'MinPeakHeight',12,'MinPeakDistance',90,'MinPeakProminence',4);
%time at peak locations
locs_t2=locs2/fs;
%Plotting after second stage of preprocessing and bandpass filtering
%Plotting detected peaks
figure();
plot(t,filPPGav2);
hold on;
plot(locs_t2,filPPGav2(locs2),'ro','MarkerFaceColor','red','MarkerSize',3);
grid on;
title({'Filtered PPG signal','Subject 1'});
xlabel('Time (s)'); 
ylabel('PPG signal[-]');


%Method 2 for actual use,
%Reference Thinh Nguyen https://uk.mathworks.com/matlabcentral/fileexchange/53364-heart-rate-spo2-using-ppg
%This one gave a slightly higher error
%It has a slightly diffferent approach using BPM estimate at each interval
%and iteratively averaging it
%compared to calculating single BPM estimate form average of all intervals
%statistical difference, E[f(x)] vs f([E(x)])
beat_cont = 0;%BPM contribution (sum) from all detected beats
beat_det = 1;%count of detected beats

%Stnadardized
%[peaks,locs]=findpeaks(filPPGav1,'MinPeakHeight',0.1,'MinPeakProminence',0.03);
for i=1:length(locs)-1
  %loops through all peaks
  %condition same as using min interval distance, reinforces
    if locs(i)<(locs(i+1)-30)
        %summing BPM estimate at interval from current peak to next
        beat_cont=beat_cont+(fs/((locs(i+1)-locs(i))))*60;
        %counting number of beats looped through
        beat_det=beat_det+1;
    end
end
if beat_det>1
    %calculating average BPM estimate by dividing summed BPM by number of
    %intervals, one less than number of beats
    avg_BPM_But_stdz= beat_cont/(beat_det-1);
end

%Unstandardized
%[peaks2,locs2]=findpeaks(filPPGav2,'MinPeakHeight',12,'MinPeakProminence',4);
for i=1:length(locs2)-1
    if locs2(i)<(locs2(i+1)-30)    
        beat_cont=beat_cont+(fs/((locs2(i+1)-locs2(i))))*60;
        beat_det=beat_det+1;
    end
end
if beat_det>1
    avg_BPM_But= beat_cont/(beat_det-1);
end

%writing to output
PPGsigOut.avg_BPM_But_stdz=avg_BPM_But_stdz;
PPGsigOut.avg_BPM_But=avg_BPM_But;

%% Part 5 Adaptive Filtering %%
% Overall parameters
FFTres = 1200;   % FFT resolution - number of bins used
WFlength = 16; % Wiener filter length

% Performance parameters
BPM_est_norm=[]; %BPM at every window for normalized signal
BPM_est=[]; % BPM at every window for unnormalized signal

% Window and iteration parameters
window   = 3 * fs;  % window length is 8 seconds
step     = 0.75 * fs;  % step size is 2 seconds
wLength = floor((length(filPPGav1)-window)/step) + 1;  % total number of windows(estimates)
    
%parameter for nested for loop for window selection
j=1;
k=1;
    for i =  [1 :  wLength]
        %selects the current segment of data based on window number
        curSeg = (i-1)*step+1 : (i-1)*step+window;
        %filtered PPG signals in our segment
        PPG_ave = filPPGav1(:,curSeg); PPG_ave2 = filPPGav2(:,curSeg);
        %fitered accelerometer signals in our segment
        ACC_X_seg = ACC_X(:,curSeg); ACC_Y_seg = ACC_Y(:,curSeg); ACC_Z_seg = ACC_Z(:,curSeg);
        
        
        % Frequency domain conversion of signals
        PPG_ave_FFT = fft(PPG_ave,FFTres);
        PPG_ave2_FFT = fft(PPG_ave2,FFTres);
        
        %acceleration FFT
        ACC_X_FFT= fft(ACC_X_seg,FFTres);
        ACC_Y_FFT= fft(ACC_Y_seg,FFTres);
        ACC_Z_FFT= fft(ACC_Z_seg,FFTres);
        
        %Reference A Temko "Accurate Heart Rate Monitoring During Physical Exercises Using PPG," 
        %in IEEE Transactions on Biomedical Engineering,vol. 64,no. 9,pp. 2016-2024,Sept. 2017
        % Wiener filtering PPG-ACC
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
        
        %calculating a window's clean signal in time-domain
        PPG_ave_IFFT_Cln(i,:)= ifft(PPG_ave_FFT.*(WF11));
        PPG_ave2_IFFT_Cln(i,:) = ifft(PPG_ave2_FFT.*(WF12));
        
        %calculating a BPM estimate for every shifted window
        %Peak detection
        [peaks3,locs3]=findpeaks(PPG_ave_IFFT_Cln(i,:),'MinPeakHeight',0.2,'MinPeakProminence',0.15);
        [peaks4,locs4]=findpeaks(PPG_ave2_IFFT_Cln(i,:),'MinPeakHeight',60,'MinPeakProminence',40);
        
        
        for k=1:length(locs3)-1
            if locs3(k)<(locs3(k+1)-90)    
            beat_cont=beat_cont+(fs/((locs3(k+1)-locs3(k))))*60;
            beat_det=beat_det+1;
            end
        end
        if beat_det>1
        BPM_est_norm(i)= beat_cont/(beat_det-1);
        end
        
        
        for k=1:length(locs4)-1
            if locs4(k)<(locs4(k+1)-90)    
            beat_cont=beat_cont+(fs/((locs4(k+1)-locs4(k))))*60;
            beat_det=beat_det+1;
            end
        end
        if beat_det>1
        BPM_est(i)= beat_cont/(beat_det-1);
        end
        
        %reconstructing the output signal with all non-overlapping windows
        %These differ by a multiple of four because of the step used
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


%error calculation for the entrie signal
[peaks5,locs5]=findpeaks(PPG_Clean_norm,'MinPeakHeight',0.15,'MinPeakProminence',0.1);
[peaks6,locs6]=findpeaks(PPG_Clean ,'MinPeakHeight',60,'MinPeakProminence',40);

%time at peak locations
locs_t5=locs5/fs;
%time at peak locations
locs_t6=locs6/fs;
        
%calculating BPM errors in the final signal normalized
        for k=1:length(locs5)-1
            if locs5(k)<(locs5(k+1)-90)    
            beat_cont=beat_cont+(fs/((locs5(k+1)-locs5(k))))*60;
            beat_det=beat_det+1;
            end
        end
        if beat_det>1
        avgBPM_SigNorm= beat_cont/(beat_det-1);
        end
        
%calculating BPM errors in the final signal       
        for k=1:length(locs6)-1
            if locs6(k)<(locs6(k+1)-90)    
            beat_cont=beat_cont+(fs/((locs6(k+1)-locs6(k))))*60;
            beat_det=beat_det+1;
            end
        end
        if beat_det>1
        avgBPM_Sig= beat_cont/(beat_det-1);
        end


%writing to output
PPGsigOut.PPG_Clean_norm=PPG_Clean_norm;
PPGsigOut.PPG_Clean=PPG_Clean;
PPGsigOut.BPM_est_norm=BPM_est_norm;
PPGsigOut.BPM_est=BPM_est;
PPGsigOut.avgBPM_SigNorm=avgBPM_SigNorm;
PPGsigOut.avgBPM_Sig=avgBPM_Sig;

t2=0:1/fs:(length(PPG_Clean)/fs);%for plotting, time values from sample rate and index
t2=t2(2:length(t2));
%Plotting figures of clean PPG
figure();
plot(t2,PPG_Clean_norm);
hold on;
plot(locs_t5,PPG_Clean_norm(locs5),'ro','MarkerFaceColor','red','MarkerSize',3);
grid on;
title({'Standardized processed PPG signal','Subject 1'})
xlabel('Time (s)');
ylabel('PPG signal [-]');
figure();
plot(t2,PPG_Clean);
hold on;
plot(locs_t6,PPG_Clean(locs6),'ro','MarkerFaceColor','red','MarkerSize',3);
grid on;
title({'Processed PPG signal','Subject 1'})
xlabel('Time (s)');
ylabel('PPG signal [-]');
xlim([5 45]);

%Plotting estimated BPM against HR ground truth normalized
figure();
plot(BPM_est_norm);
hold on;
plot(BPM_est);
grid on;
title({'BPM estimate comparison on standardization and normalization of signal','Subject 1'})
xlabel('Window number');
ylabel('BPM estimate');
legend('Estimated BPM (norm. sig)','Estimated BPM');



