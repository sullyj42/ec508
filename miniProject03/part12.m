close all; clear all
SNRdb = 0; SNR = 10^(SNRdb/10); 
N = 1; P = SNR*N;  
% Let's generate 100 channel realizations, 100 times
% Given h is iid gaus(0,1)
M = 1E5; L = 1; 
h       = 1/sqrt(2) * randn(L,M) + 1j/sqrt(2) * randn(L,M);
[ht, I] = sort(abs(h).^2, 2); 
hd      = mean(ht,1); 
effNL      = sort(N./hd);
stem(effNL(1:00)); xlabel('h[m]'); ylabel('N/h[m]'); set(gca, 'YScale', 'log')
effPower = 0; i = 1; 
while effPower < P*M && i < length(effNL)
    effPower = effPower + (effNL(i+1) - effNL(i))*i; % This is some janky *ss code. 
    wl = effNL(i+1);                              % This is not exact, but who cares? 
    i = i+1; 
end
if i == length(h) % Did not fill all the buckets
   warning('Yo this aint tested')
   wl = wl + (P-effPower)./length(h);
else 
    lastLevel = i-1;
    wl = effNL(lastLevel); % Take a step back so you are below the power constraint
                  % This could be made exact, but... this is hard
                  % Let's try for lolz... 
    effPowerOld  = effPower -  (effNL(i) - effNL(i-1)).*(i-1); 
    partialLevel = (P*M - effPowerOld)./(i-1);
    wl           = effNL(i-1) + partialLevel; 
    err = sum(wl - effNL(1:i-1)) -  P*M;
    if 1E-8 < err
        warning(['Power level is off by ', num2str(err)])
    end
end


clf; stem(effNL(1:100)); xlabel('h[m]'); ylabel('N/h[m]'); set(gca, 'YScale', 'log'); 
hold on; plot([0, 100], [wl, wl], 'o-'); hold off;
legend('Average Effective Noise Levels','Average Effective Water Level')

wfUnsort = 1./abs(h(1,:)).^2; % Take one random slice
stem(wfUnsort(1:100)); xlabel('h[m]'); ylabel('N/h[m]'); set(gca, 'YScale', 'log'); 
hold on; plot([0, 100], [wl, wl], 'o-'); hold off; 
legend('Instantaneous Effective Noise Levels', 'Average Effective Water Level');
figure; plot(effNL, '.-'); hold on; plot([0, length(effNL)], [wl, wl], 'o-'); set(gca, 'YScale', 'log')
xlabel('Channel Instance'); ylabel('Effective Noise Levels'); title('Overall Simulation'); 
xlim([0, i+i*0.1]); 
%% Part B

SNRdb = -20:20; SNR = 10.^(SNRdb/10); 
N = 1;  
% Let's generate 100 channel realizations, 100 times
% Given h is iid gaus(0,1)
M = 1E5; L = 1; 
h       = 1/sqrt(2) * randn(L,M) + 1j/sqrt(2) * randn(L,M);
[ht, I] = sort(abs(h).^2, 2); 
hd      = mean(ht,1); 
effNL      = sort(N./hd);

for ii = 1:length(SNR)
 P = SNR(ii)*N; 
effPower = 0; i = 1; 
while effPower < P*M && i < length(effNL)
    effPower = effPower + (effNL(i+1) - effNL(i))*i; % This is some janky *ss code. 
    wl = effNL(i+1);                              % This is not exact, but who cares? 
    i = i+1; 
end
if i == length(h) % Did not fill all the buckets
   warning('Yo this aint tested')
   wl = wl + (P-effPower)./length(h);
else 
        lastLevel = i-1;
        wl = effNL(lastLevel); % Take a step back so you are below the power constraint
                  % This could be made exact, but... this is hard
                  % Let's try for lolz... 
        effPowerOld  = effPower -  (effNL(i) - effNL(i-1)).*(i-1); 
        partialLevel = (P*M - effPowerOld)./(i-1);
        wl           = effNL(i-1) + partialLevel; 
        err = sum(wl - effNL(1:i-1)) -  P*M;
        if 1E-8 < err
            warning(['Power level is off by ', num2str(err)])
        end
end
Pm = (wl-effNL(1:(i-1))); 
h = sort(h); 
C1(ii) = sum(log2(1 + Pm.*abs(h((end-length(Pm)+1):end)).^2./N))./M;  % Can't use mean cuz average over M...
C2(ii) = mean(log2(1+(abs(h)).^2.*SNR(ii)));
end
awgn = log2(1 + std(h) .* SNR./N); 
hold on; plot([0, 100], [wl, wl], 'o-'); hold off;
legend('Average Effective Noise Levels','Average Effective Water Level')

figure; 
plot(SNRdb, awgn); hold on; 
plot(SNRdb, C1)
plot(SNRdb, C2) 
xlabel('SNR (dB)'); ylabel('Capacity (bits/s/hz)'); 
legend('AWGN', 'Full CSI', 'CSIR')
%% Part C
figure; xlabel('SNR (dB)');  ylabel('C/C_{awgn}') 
plot(SNRdb, C1./awgn); hold on; 
plot(SNRdb, C2./awgn); 
xlim([-20, 10]); legend('Full CSI', 'CSIR'); 


