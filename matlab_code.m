clear all;
clc;                                                                        %house keeping
load mit200.mat                                                             % ECG data
wt = modwt(ecgsig,5);                                                       % maximal overlap discrete wavelet transform
wtrec = zeros(size(wt));
wtrec(4:5,:) = wt(4:5,:);                                                   
y = imodwt(wtrec,'sym4');                                                   %reconstruction signal usng Inverse maximal overlap discrete wavelet transform
fs=360;
[qrspeaks,locs] = findpeaks(y,tm,'MinPeakHeight',0.35,...
    'MinPeakDistance',0.150);                                               % finding r peaks
%[qspeaks,loc_qs] = findpeaks(-y,tm,'MinPeakHeight',0.14,'MinPeakDistance',0.15); %% finding q and s peaks
[qspeaks,loc_qs] = findpeaks(-y,tm,'MinPeakHeight',0.14,'MinPeakProminence',0.2); %% finding q and s peaks
[speaks,loc_s] = findpeaks(-y,tm,'MinPeakHeight',0.35,'MinPeakProminence',0.15); %% finding s peaks

j=1;
for i=1:length(loc_qs)                                                      % separating q peaks from s peaks
    if -qspeaks(i)>-0.35
        qpeaks(j)=qspeaks(i);
        loc_q(j)=loc_qs(i);
        j=j+1;
    end
end

figure;                                                                     %graphs for wavelet based detection
plot(tm,y)
hold on
plot(loc_q,-qpeaks,'ko')
hold on;
plot(locs,qrspeaks,'r*')
hold on
plot(loc_s,-speaks,'b*')
legend('ECG','Q peaks','R peaks', 'S peaks')
xlabel('Seconds')
title('QRS Peaks Localized by Wavelet Transform')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);

figure;
plot(tm,ecgsig);
hold on;
plot(ann./fs,ecgsig(ann),'r*');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);

%%  pan tompkins
%%bandpass filter
% [ecg_bp d]=bandpass(ecgsig,[5,15],250);
Wn=[5 15]*2/fs;                                                           % cutt off based on fs
 N = 3;                                                                     % order of 3 less processing
[a,b] = butter(N,Wn);                                                      % bandpass filtering
ecg_bp = filtfilt(a,b,ecgsig);
ecg_bp = ecg_bp/ max( abs(ecg_bp));                                         %normalization
figure;
plot(ecg_bp);
hold on;
plot(ecgsig,'r');
legend('Filtered ECG','Unfiltered ECG')
xlabel('Seconds')
title('Bandpass filtered ECG over Unfiltered ECG')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);

%% derivative filter
 num = [1 2 0 -2 -1].*(1/8)*fs;  
 den=1;
 ecg_d = filtfilt(num,1,ecg_bp);
 d_fil = ecg_d/max(ecg_d);

figure;
plot(tm,ecg_bp);
hold on;
plot(tm,d_fil,'r');
xlabel('Seconds')
legend('Bandpass Filtered ECG','Derivative filtered ECG')
title('Derivative Filter')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);

%%
 %freqz(num,den,10000,360)

%%
delay=0;
ecg_sq=(d_fil).^2;
figure;
plot(tm,ecg_sq);
xlabel('Seconds')
title('Squared Signal')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);

%% Moving average filter
% ecg_m = conv(ecg_sq ,ones(1 ,round(0.150*360))/round(0.150*360));
win=(1/round(0.150*360))*ones(1 ,round(0.150*360));                         %number of samples to average is chosen in order to average on windows of 150 ms                 
ecg_m=filter(win,1,ecg_sq);
delay = delay + round(0.150*360)/2;
%ecg_ma=filter(37.5,1,ecg_sq);

figure;
plot(ecg_m);
legend('Moving Average Filtered signal')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);

%%
figure;
plot(ecg_bp);
hold on;
plot(ecg_m);
legend('Bandpass Filtered signal','Moving Average Filtered signal')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);


%%
tm = tm(1:end-delay);
ecg_bp = ecg_bp(1:end-delay);                                              % Removing delay due to MA filter

sf = ecg_m;
sf(1:delay) = [];
ecg_m=sf;
%%
figure;                                                                     %figures without delay
plot(ecg_bp);
hold on;
plot(ecg_m);
legend('Bandpass Filtered signal','Moving Average Filtered signal')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);


%% Thresholding
[pks,locs] = findpeaks(ecg_m,'MINPEAKDISTANCE',round(0.2*360));             % finding peaks of integrated signal
[pks_f,locs_f] = findpeaks(ecg_bp,'MinPeakHeight',0.14,'MINPEAKDISTANCE',round(0.2*360)); % finding peaks of bandpassed filtered signal
figure;
plot(ecg_bp)
hold on;
plot(locs_f,pks_f,'r*');
hold on;
plot(ecg_m)
hold on;
plot(locs,pks,'k*');
legend('Bandpass Filtered signal','Moving Average Filtered signal','R peak BP signal','Peak MA Signal')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);




%%
LLp = length(pks);
THR_SIG = max(ecg_m(1:2*fs))*1/3;                                          % 0.25 of the max amplitude 
THR_NOISE = mean(ecg_m(1:2*fs))*1/2;                                       % 0.5 of the mean signal is considered to be noise
SIG_LEV= THR_SIG;
NOISE_LEV = THR_NOISE;
thr=NOISE_LEV+0.25*(SIG_LEV-NOISE_LEV);

THR_SIG1 = max(ecg_bp(1:2*fs))*1/3;                                          % 0.25 of the max amplitude 
THR_NOISE1 = mean(ecg_bp(1:2*fs))*1/2; 
SIG_LEV1 = THR_SIG1;                                                        % Signal level in Bandpassed filter
NOISE_LEV1 = THR_NOISE1;
THR_F=NOISE_LEV1+0.25*(SIG_LEV1-NOISE_LEV1);



j=1;
k=1;
for i=1:LLp
    if pks(i)>=thr
        SIG_LEV=0.125*pks(i)+0.875*SIG_LEV;
        if pks_f(k)>=THR_F
           SIG_LEV1=0.125*pks_f(k)+0.875*SIG_LEV1;
           qrs_peak(j)=pks(i);                                              %% implement another threshold conditon using bandpass signal
           qrs_locs(j)=locs(i);
           j=j+1;
           k=k+1;
        else
             NOISE_LEV1=0.125*pks_f(i)+0.875*NOISE_LEV1;
        end
        
    else
        NOISE_LEV=0.125*pks(i)+0.875*NOISE_LEV;
    end
    THR(i)=thr;
    thr=NOISE_LEV+0.25*(SIG_LEV-NOISE_LEV);                                 %Updating the threshold
    THR_F=NOISE_LEV1+0.25*(SIG_LEV1-NOISE_LEV1);

    
end
%%
figure;
plot(ecg_m);
hold on;
plot(locs,THR);
hold on;
plot(qrs_locs,qrs_peak,'r*');
legend('Moving Average Filtered Signal','Threshold','QRS Complex')
title('QRS Peaks Localized by Pan-Tompkins');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
%%
figure;
plot(ecg_m);
hold on;
plot(qrs_locs,qrs_peak,'k*');
hold on;
plot(ann,ecg_m(ann),'r*');
legend('Integrated Signal','Pan-Tompkins Found Ann','Origional Ann')
title('QRS Peaks Localized by Pan-Tompkins');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);

%%
figure;
plot(ecgsig);
hold on;
plot(qrs_locs,qrs_peak,'k*');
hold on;
plot(ann,ecgsig(ann),'r*');
legend('Origional ECG','Pan-Tompkins Found Ann','Origional Ann')
title('QRS Peaks Localized by Pan-Tompkins');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);

%%


[pks_raw,locs_raw] = findpeaks(ecg_bp,'MINPEAKDISTANCE',round(0.2*360));

figure;
plot(ecg_bp);
hold on;
plot(locs_raw,pks_raw,'r*');
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);


%%
qrs_n=length(qrs_locs);
num_beats=length(qrs_locs);

HR = num_beats*60/27.775; %% heart rate based on Pan-Tompkins
