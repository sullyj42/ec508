clc; close all; clear all; 
% Sections 3.1, 3.2, and 3.3 of Tse & Viswanat
%% Error probability sending BPSK symbols
%{

SNRdb   = linspace(0,10, 10);
SNR     = 10.^(SNRdb/10); L = 1:4;
L       = 1:3; 
N       = 1000; 
err     = zeros(length(L), length(SNR), N); 
errSum  =  zeros(length(L), length(SNR)); 
%{
for i = 1:N
    err(:,:, i) = calcErr(SNR, L);
end
while ~isequal(zeros(length(L), length(SNR)), (errSum == 0))
    i = i+1; err(:,:,i) = calcErr(SNR, L); errSum = mean(err, 3);
end
%}

%% Part B: 

[P, Pe] = calcExact(SNR, L); 

for iL = 1:length(L)
    semilogy(SNRdb, (errSum(iL, :)), '.-k'); hold on; 
    semilogy(SNRdb, (P(iL, :))  , '.-r'); 
    semilogy(SNRdb, (Pe(iL, :)) , '.-b'); 
end
    legend('Simulation', 'Exact', 'Estimate'); 
    xlabel('SNR (dB)'); ylabel('Probability of Error'); 
    title('BPSK Error Rate for L = 1:3'); 
    print('partB', '-dpng'); 
%}    
%% Part C:
    clearvars; close all; 
    SNRdb   = 0;
    SNR     = 10.^(SNRdb/10); 
    L       = 1:3; 
    N       = 10000; 
    M       = 20; Mc = 100;  mkdir(sprintf('N_%01.0f_Mc_%01.0f/', N, Mc));
    err     = zeros(length(L), length(SNR), N, M); 
    errSum  = zeros(length(L), length(SNR)); 
    errESum = errSum; 
    [P, Pe] = calcExact(SNR, L); 

for m= 1:M % Iterate over the number of training symbols 
    for i = 1:N % Iterate over the number of trials
        [err(:,:, i, m), errE(:,:,i,m)] = calcErr_Hest(SNR, L, m, Mc) ; % Store the array for each training symbol 
    end
    errSum(:, :, m)  = mean(err(:,:,:,m), 3);
    errESum(:, :, m) = mean(errE(:,:,:,m), 3); 
end
figure; 
if length(SNR)==1 
    for i = 1:length(L)
        a = reshape(errSum(i, 1, :), 1, M); 
        b = reshape(errESum(i, 1, :), 1, M); 
        semilogy(1:M, a, '.-k'); hold on; 
        semilogy(1:M, b, '.-b'); 
        semilogy([1, M], [P(i), P(i)], 'o-r'); 
    end
    legend('Measured Channel', 'Exact Channel', 'Theoretical'); 
    xlabel('M_{train}'); ylabel('P_{error}'); 
    title(sprintf('Probability of Error for various channel knowledge\nSNR: %01.0fdB, N: %01.0f, Mc %01.0f', [0 N, Mc])); 
    print(sprintf('SNR%01.0f_N%01.0f_Mc%01.0f', [SNR N, Mc]), '-dpng')
end

%% Part D:
    clearvars; close all; 
    SNRdb   = 0:10;
    SNR     = 10.^(SNRdb/10); 
    L       = 1:3; 
    N       = 2000; 
    M       = 5; Mc = 100;  mkdir(sprintf('N_%01.0f_Mc_%01.0f/', N, Mc));
    err     = zeros(length(L), length(SNR), N, M); 
    errSum  = zeros(length(L), length(SNR)); 
    errESum = errSum; 
    [P, Pe] = calcExact(SNR, L); 

for m= 1:M % Iterate over the number of training symbols 
    parfor i = 1:N % Iterate over the number of trials. Not sure if parfor actually does anything here
        [err(:,:, i, m), errE(:,:,i,m)] = calcErr_Hest(SNR, L, m, Mc) ; % Store the array for each training symbol 
    end
    errSum(:, :, m)  = mean(err(:,:,:,m), 3);
    errESum(:, :, m) = mean(errE(:,:,:,m), 3);
    
    % Plot for each diversity level
    for iL = 1:length(L)
        semilogy(SNR, (errSum(iL, :, m)), '.-k'); hold on; 
        semilogy(SNR, (errESum(iL, :, m)), '.-b'); hold on; 
        semilogy(SNR, (P(iL, :)), '.-r'); 
    end   
    
        legend('Simulation (measured channel)', 'Simulation (exact channel)', 'Exact (theoretical)'); 
        xlabel('SNR (dB)'); ylabel('Probability of Error'); 
        title(['BPSK Error Rate for estimated Channel, M_{train} = ', num2str(m)]); 
        hold off; drawnow; pause(0.5)
        print(sprintf('N_%01.0f_Mc_%01.0f/m', N, Mc, m), '-dpng');
    %}
end

    
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
        y   = h * x  + w;                 % Scale by the channel effects and add the random noise
%{
         yh  = ((y.*hh)); 
         yhn = sum(yh); % This seems wrong
         if      yhn >= 0;    xh = 1; 
         elseif  yhn <  0;    xh = -1; 
         else;               warning('bad'); 
         end  
%}
        if norm(y - h*sqrt(SNR)) <= norm(y + h*sqrt(SNR)) 
            xh = 1; 
        elseif norm(y - h*sqrt(SNR)) >= norm(y + h*sqrt(SNR))
            xh = -1;       
        else; warning("You're dumb"); 
        end % Perform ML Assignment
        
        err(iL, iSNR) = ~isequal(xh, a);
    end
end
end

function [err, errE] = calcErr_Hest(SNR, L, M, Mc) 
    % Determines symbol based on measured (err) vs. exact (errE) channel knowledge
    err  = zeros(length(L), length(SNR));
    errE = err;  
    for iL = 1:length(L)

        % Generate the true channel effect
        h   = 1/sqrt(2) * randn(L(iL),1) + 1i/sqrt(2) * randn(L(iL),1); 
        hh  = conj(h)./(abs(h));
        for iSNR = 1:length(SNR)


            % "Preamble / Pilot"
            A = sqrt(SNR(iSNR))/(1 + M * SNR(iSNR)); 
            % clearvars B % Varying length array, do not want to 
            B = zeros(length(h), length(M));  
            % Measure/estimate the channel effects
            for m = 1:M
                w    = 1/sqrt(2) * randn(L(iL),1) + 1i/sqrt(2) * randn(L(iL),1); 
                B(:, m) = h.*sqrt(SNR(iSNR)) + w;
            end
            hg = A.*sum(B,2);

            for mm = (M+1):(Mc)

                w   = 1/sqrt(2) * randn(L(iL),1) + 1i/sqrt(2) * randn(L(iL),1); 
                a   = (2*randi([0 1]) - 1);       % Produces +/-
                x   = sqrt(SNR(iSNR))*a;          % Scale the symbol by SNR   --> This could be vectorized

                y   = h * x  + w;                 % Scale by the channel effects and add the random noise
                % These evaluations are suppppppper slow. 
                % Taking 60% of time (profiler)
                if norm(y - hg*sqrt(SNR), 1) <= norm(y + hg*sqrt(SNR),1) % Perform ML Assignment for measured channel
                    xh = 1; 
                elseif norm(y - hg*sqrt(SNR),1 ) >= norm(y + hg*sqrt(SNR),1)
                    xh = -1;       
                else; warning("You're dumb"); 
                end 

                % These evaluations are suppppppper slow. 
                % Taking 60% of time (profiler)
                if norm(y - h*sqrt(SNR),1) <= norm(y + h*sqrt(SNR),1) % Perform ML Assignment for exact channel
                    xhh = 1; 
                elseif norm(y - h*sqrt(SNR),1) >= norm(y + h*sqrt(SNR),1)
                    xhh = -1;       
                else; warning("You're dumb"); 
                end 

                err(iL, iSNR)  = err(iL, iSNR)  + ~isequal(xh, a);
                errE(iL, iSNR) = errE(iL, iSNR) + ~isequal(xhh, a);

            end
            err(iL, iSNR)  = err(iL, iSNR) ./((Mc-M)); % Find the average error for those sections
            errE(iL, iSNR) = errE(iL, iSNR)./((Mc-M)); 
            % In each Mc-interval, we evaluated Mc-M symbols
        end
    end
end

function [P, Pe] = calcExact(SNR, L) 
    % Constructs two LxSNR matrices 
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
end

