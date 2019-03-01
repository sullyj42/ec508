clc; close all; clear all; 
% Sections 3.1, 3.2, and 3.3 of Tse & Viswanat
%% Error probability sending BPSK symbols


SNRdb   = linspace(0,10, 10);
SNR     = 10.^(SNRdb/10); L = 1:4;
L       = 1:3; 
N       = 100000; 
err     = zeros(length(L), length(SNR), N); 
errSum  =  zeros(length(L), length(SNR)); 

for i = 1:N
    err(:,:, i) = calcErr(SNR, L);
end
while ~isequal(zeros(length(L), length(SNR)), (errSum == 0))
    i = i+1; err(:,:,i) = calcErr(SNR, L); errSum = mean(err, 3);
end

%% Part B: 

for iL = 1:length(L)
    l = 0:(L(iL)-1); 
    for iSNR = 1:length(SNR)
        u = sqrt(SNR(iSNR)/(1+SNR(iSNR))); 
        for l = 0:(L(iL) - 1)
            temp(l+1)  = nchoosek(L(iL)-1+l, l) * ((1 + u)/2)^l;
        end
        P(iL, iSNR)  = sum(temp)*((1-u)/2)^L(iL); 
        Pe(iL, iSNR) = nchoosek(2*L(iL)-1, L(iL))/(4*SNR(iSNR))^L(iL); 
    end
end

for iL = 1:length(L); 
    semilogy(SNRdb, (errSum(iL, :)), '.-k'); hold on; 
    semilogy(SNRdb, (P(iL, :))  , '.-r'); 
    semilogy(SNRdb, (Pe(iL, :)) , '.-b'); 
end
    legend('Simulation', 'Exact', 'Estimate'); 
    xlabel('SNR (dB)'); ylabel('Probability of Error'); 
    title('BPSK Error Rate for L = 1:3'); 
%% Supporting functions

function err = calcErr(SNR, L) 
err = zeros(length(L), length(SNR)); 
for iL = 1:length(L)
    h   = 1/sqrt(2) * randn(L(iL),1) + 1i/sqrt(2) * randn(L(iL),1); 
    hh  = conj(h)./(abs(h));
    for iSNR = 1:length(SNR)
        w   = 1/sqrt(2) * randn(L(iL),1) + 1i/sqrt(2) * randn(L(iL),1); 
        a   = (2*randi([0 1]) - 1);       % Produces +/-
        x   = sqrt(SNR(iSNR))*a;          % Scale the symbol by SNR   --> This could be vectorized
        y   = h * x  + w;              % Scale by the channel effects and add the random noise
%         yh  = ((y.*hh)); 
%         yhn = sum(yh); % This seems wrong
%         if      yhn >= 0;    xh = 1; 
%         elseif  yhn <  0;    xh = -1; 
%         else;               warning('bad'); 
%         end        
        if norm(y - h*sqrt(SNR)) <= norm(y + h*sqrt(SNR)) 
            xh = 1; 
        elseif norm(y - h*sqrt(SNR)) >= norm(y + h*sqrt(SNR))
            xh = -1; 
        else 
            warning("You're dumb")
        end
        err(iL, iSNR) = ~isequal(xh, a);
    end
end
end


