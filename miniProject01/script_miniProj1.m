% As with all good things in engineering... best intentions to take an
% object-oriented approach to keep things neat&tidy... 
% But naturally copy and paste wins out in the end
clear all; close all; 
%% Initial setup
%%%%%%%%% Cartesian Coordinates %%%%%%%%%
%  = [x0, y0; Position t(1) // refl(1)  % 
%     x1, y1; Position t(2) // refl(2)  %
%      .   .                            %
%      .   .                            %
%      .   .                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tx.pos = [0, 0]; % Transmitter fixed at origin
% Parameters that will not change
global C;       C       = 3E8;     
global fc;      fc      = 900E6; 
global lambdaC; lambdaC = C/fc; 
global Grx;     Grx     = 1; 
global Gtx;     Gtx     = 1; 
global W;       W       = 1E7;
                       % I hate globals, but this is so much easier...
load('mp1data.mat');   % Reflector Positions

% Reorganize because I am so particular
reflectors{1}.pos = dreflect1; 
reflectors{2}.pos = dreflect2; 

%% (a): Geometric set-up
% Assume a given initial reciever position and speed
Rx.posInit = [0, 0]; % Start colocated and and then spread
Rx.V       = [1, 1]; % at 1 m/s in both directions
timeArray  = linspace(0, 3600/2); % Take a half-hour of data
% Not exactly a function of time, but... kind-of
    Rx.pos     = recieverPosition(Rx.posInit, Rx.V, timeArray); 
    t_ind = length(Rx.pos); % let's look torwards end of data
    reflectors{1}.range  = pathLength(Tx.pos, Rx.pos(t_ind,:), reflectors{1}.pos);
    reflectors{2}.range  = pathLength(Tx.pos, Rx.pos(t_ind,:), reflectors{2}.pos);


% This is ugly... not sure why plot doesn't work with a [n x 2] array... 
    plot(Rx.pos(:,1), Rx.pos(:,2), '.'); hold on; plot(Tx.pos(1), Tx.pos(2), 'o'); 
    subplot(2,2,1);  title('Reflector Set 1'); 
    plot(reflectors{1}.pos(:,1), reflectors{1}.pos(:,2), '.'); hold on; 
    plot(Rx.pos(:,1), Rx.pos(:,2), '.-'); 
    xlabel('X-Coordinate'); ylabel('Y-Coordinate'); 
    title('Scenario Geometry'); legend('Reflectors','Reciever')

    subplot(2,2,2); title('Reflector Set 2'); 
    plot(reflectors{2}.pos(:,1), reflectors{2}.pos(:,2), '.'); hold on; 
    plot(Rx.pos(:,1), Rx.pos(:,2), '.-'); 
    xlabel('X-Coordinate'); ylabel('Y-Coordinate'); 
    title('Scenario Geometry'); legend('Reflectors', 'Reciever')

    subplot(2,2,3); plot(reflectors{1}.range); title('Path length to reciever at end (DP ~2500'); 
    xlabel('Reflector Number (...meaningless?)'); ylabel('Path length'); 

    subplot(2,2,4); plot(reflectors{2}.range); title('Path length to reciever at end (DP ~ 2500') 
    xlabel('Reflector Number (...meaningless?)'); ylabel('Path length'); 
% Looks as expected...
    
%% (b): calculate, at a given time t, the attenuation a_i(t) and delay T_i(t) for each path
% Requires pathlength (range)
clear Rx.pos timeArray; 
Rx.posInit = [100, 100]; % Arbitray;
Rx.V       = [1, 1]; 
timeArray   =  1; % Let's take one time value; 
% May throw this into a function for clarity as a function of time; 
    Rx.pos     = recieverPosition(Rx.posInit, Rx.V, timeArray); 
    reflectors{1}.range  = pathLength(Tx.pos, Rx.pos, reflectors{1}.pos);
    reflectors{1}.rangeS = sort(pathLength(Tx.pos, Rx.pos, reflectors{1}.pos));
    reflectors{2}.range  = pathLength(Tx.pos, Rx.pos, reflectors{2}.pos); 
    reflectors{2}.rangeS = sort(reflectors{2}.range); % For neater plotting for large N reflectors

reflectors{1}.atten  = pathAttenuation(Grx, Gtx, lambdaC, reflectors{1}.range);
reflectors{1}.attenS = pathAttenuation(Grx, Gtx, lambdaC, reflectors{1}.rangeS);
reflectors{2}.atten  = pathAttenuation(Grx, Gtx, lambdaC, reflectors{2}.range);
reflectors{2}.attenS = pathAttenuation(Grx, Gtx, lambdaC, reflectors{2}.rangeS);
figure; subplot(2,1,1); plot(reflectors{1}.rangeS, reflectors{1}.attenS, '.-'); 
xlabel('Reflector Path Length'); ylabel('Attenuation (linear, A/A unitless)'); 
title('Part b: Path loss from reflectors (set 1)'); 
subplot(2,1,2);  plot(reflectors{2}.rangeS, reflectors{2}.attenS, '.-');
xlabel('Reflector Path Length'); ylabel('Attenuation (linear, A/A unitless)'); 
title('Reflectors (set 2)'); 

%% (c): Continuous time channel IR and frequency response
% Plot power spectrum in decibels at a given time
% I do not understand this part and my results are garbage... moving on

% I am really not sure how to define these variables as a function of time
% I think I need to get more analytical

reflectors{1}.pDelay = pathDelay(reflectors{1}, C); 
reflectors{2}.pDelay = pathDelay(reflectors{2}, C);
% It seems silly to plot anything besides 
[reflectors{1}.h, channelImpt] = channelImpResp(reflectors{1}, reflectors{1}.pDelay, length(reflectors{2}.pDelay)+100); 
figure; stem(channelImpt, reflectors{1}.h); figure; 
fs = 1/(channelImpt(2)-channelImpt(1)); f = fs*(-length(channelImpt)/2:length(channelImpt)/2)./length(channelImpt); 
Y = fft(reflectors{1}.h); figure; plot(f(1:end-1), abs(fftshift(Y))); title('Impulse response for IR of reflector set 1'); 

[reflectors{2}.h, channelImpt] = channelImpResp(reflectors{2}, reflectors{2}.pDelay, length(reflectors{2}.pDelay)+100); 
figure; stem(channelImpt, reflectors{2}.h); figure;
fs = 1/(channelImpt(2)-channelImpt(1)); f = fs*(-length(channelImpt)/2:length(channelImpt)/2)./length(channelImpt); 
Y = fft(reflectors{2}.h); figure; plot(f(1:end-1), abs(fftshift(Y))); title('FT response for IR of reflector set 2'); 



%% (d): Discrete time baseband equivalent model
% h = @(m, w, rx) (sum( rx.
figure; 
h = dtbIR(reflectors{1}); 
plot(abs(h(1,:)));  figure; clear h;
xlabel('Taps'); ylabel('Amplitude?'); 
h = dtbIR(reflectors{2}); 
plot(abs(h(1,:))); xlabel('Taps'); ylabel('Amplitude?'); 
%% (e): Delay spread
reflectors{1} = delaySpread(reflectors{1}); 
reflectors{2} = delaySpread(reflectors{2}); 

%% (f): Doppler spread
% This is more tricky. Have to account for v and directions
% Wouldn't this actually be a bistatic doppler calculation??

Rx.V       = [1, 1];     Rx.posInit = [0 0]; 
Rx.pos     = recieverPosition(Rx.posInit, Rx.V, 0:10);
reflectors{1} = dopplerSpread(reflectors{1}, Rx); 
reflectors{2} = dopplerSpread(reflectors{2}, Rx); 

%% (g): ...Load data. Already done
%% (h & i): Parametric evaluation
fc = 900E6; 
W  = 1E6; 
% Take one: r1
clear Rx; clear reflectors; 
reflectors.pos = dreflect1; 
Rx.V       = [1,1]; Rx.posInit = [1000, 900]; 
Rx.pos     = recieverPosition(Rx.posInit, Rx.V, 0); 
reflectors.range  = pathLength(Tx.pos, Rx.pos, reflectors.pos);
reflectors.atten  = pathAttenuation(Grx, Gtx, lambdaC, reflectors.range);
reflectors.pDelay = pathDelay(reflectors, C);
[reflectors.h, channelImpt] = channelImpResp(reflectors, reflectors.pDelay, length(reflectors.pDelay)+2000); 
figure; stem(channelImpt, reflectors.h); xlabel('Time (seconds'); 
fs = 1/(channelImpt(2)-channelImpt(1)); f = fs*(-length(channelImpt)/2:length(channelImpt)/2)./length(channelImpt); 
Y = fft(reflectors.h); figure; plot(f(1:end-1), abs(fftshift(Y))); title('Fourier transform before shift')
xlim([-10E6, 10E6]); xlabel('Frequency (Hz) (before shift)'); 
title('Frequency domain of impulse response (ideally flat)'); 
delaySpread(reflectors); dopplerSpread(reflectors, Rx); 
% Take two: r1
clear Rx reflector
Rx.V       = [10,10]; 
Rx.posInit = [500, 1200]; 
reflectors.pos = dreflect1; 
Rx.pos     = recieverPosition(Rx.posInit, Rx.V, 0); 
reflectors.range  = pathLength(Tx.pos, Rx.pos, reflectors.pos);
reflectors.atten  = pathAttenuation(Grx, Gtx, lambdaC, reflectors.range);
reflectors.pDelay = pathDelay(reflectors, C);
[reflectors.h, channelImpt] = channelImpResp(reflectors, reflectors.pDelay, length(reflectors.pDelay)+2000); 
figure; stem(channelImpt, reflectors.h); xlabel('Time (seconds'); 
fs = 1/(channelImpt(2)-channelImpt(1)); f = fs*(-length(channelImpt)/2:length(channelImpt)/2)./length(channelImpt); 
Y = fft(reflectors.h); figure; plot(f(1:end-1), abs(fftshift(Y))); title('Fourier transform before shift')
xlim([-10E6, 10E6]); xlabel('Frequency (Hz) (before shift)'); 
title('Frequency domain of impulse response (ideally flat)'); 
delaySpread(reflectors); dopplerSpread(reflectors, Rx); 

% Take three: r2
clear Rx reflector
Rx.V       = [1,1]; 
Rx.posInit = [1000, 900]; 
reflectors.pos = dreflect2; 
Rx.pos     = recieverPosition(Rx.posInit, Rx.V, 0); 
reflectors.range  = pathLength(Tx.pos, Rx.pos, reflectors.pos);
reflectors.atten  = pathAttenuation(Grx, Gtx, lambdaC, reflectors.range);
reflectors.pDelay = pathDelay(reflectors, C);
[reflectors.h, channelImpt] = channelImpResp(reflectors, reflectors.pDelay, length(reflectors.pDelay)+2000); 
figure; stem(channelImpt, reflectors.h); xlabel('Time (seconds'); 
fs = 1/(channelImpt(2)-channelImpt(1)); f = fs*(-length(channelImpt)/2:length(channelImpt)/2)./length(channelImpt); 
Y = fft(reflectors.h); figure; plot(f(1:end-1), abs(fftshift(Y))); title('Fourier transform before shift')
xlim([-10E6, 10E6]); xlabel('Frequency (Hz) (before shift)'); 
title('Frequency domain of impulse response (ideally flat)'); 
delaySpread(reflectors); dopplerSpread(reflectors, Rx); 

% Take four: r2
clear Rx reflector
Rx.V       = [10,10]; 
Rx.posInit = [500, 1200]; 
reflectors.pos = dreflect2; 
Rx.pos     = recieverPosition(Rx.posInit, Rx.V, 0); 
reflectors.range  = pathLength(Tx.pos, Rx.pos, reflectors.pos);
reflectors.atten  = pathAttenuation(Grx, Gtx, lambdaC, reflectors.range);
reflectors.pDelay = pathDelay(reflectors, C);
[reflectors.h, channelImpt] = channelImpResp(reflectors, reflectors.pDelay, length(reflectors.pDelay)+2000); 
figure; stem(channelImpt, reflectors.h); xlabel('Time (seconds'); 
fs = 1/(channelImpt(2)-channelImpt(1)); f = fs*(-length(channelImpt)/2:length(channelImpt)/2)./length(channelImpt); 
Y = fft(reflectors.h); figure; plot(f(1:end-1), abs(fftshift(Y))); title('Fourier transform before shift')
xlim([-10E6, 10E6]); xlabel('Frequency (Hz) (before shift)'); 
title('Frequency domain of impulse response (ideally flat)'); 
delaySpread(reflectors); dopplerSpread(reflectors, Rx); 

%% (i): 



%% functions
function range = pathLength(Tx, Rx, reflectors)
    % Assumes fixed transmitter and reciever
    % Fixed, but multiple reflectors
    % Calculates round-trip trip
    if ~isequal(reflectors(1,:), Tx) || ~isequal(reflectors(1,:), Rx) % Calculate the direct path; 
        reflectors = [0, 0; reflectors]; 
    end
    range = zeros(1, length(reflectors)); 
    
    for i = 1:length(reflectors)
        temp1    = sqrt((Tx(1) - reflectors(i, 1)).^2 + (Tx(2) - reflectors(i, 2)).^2);
        temp2    = sqrt((Rx(1) - reflectors(i, 1)).^2 + (Rx(2) - reflectors(i, 2)).^2);
        range(i) = temp1 + temp2;
    end
end

function Rx = recieverPosition(RxInitial, V, t) 
    % Takes initial reciever position and estimates new position
    % Assumes V is time-invariant
    % Assumes t is an array of times
    
    Rx = zeros(size(V) + [1, 0]);
    
    Rx(1,:) = RxInitial; 
    for i = 1:length(t)
        Rx(i+1, :) = RxInitial + V.*t(i); 
    end
end

function pathAtten   = pathAttenuation(Grx, Gtx, lambdaC, pathLength)
    % Strictly speaking, I would not consider Tx/Rx gains
    % attenuation
    % Pathlength is allowed to change with time
    % All other variables should be time invariant, scalars
    pathAtten = ( sqrt(Gtx*Grx)*lambdaC/(4*pi) ) ./ (pathLength); 
end

function pDelay = pathDelay(rx, C) % Global variables suck
    pLength = rx.range; 
    % Returns delay time in seconds
    pDelay      = pLength ./ (C) ; 
end

function [channelImp channelImpt] = channelImpResp(rx, t, N)
    % Expects rx to be a struct with precalculated delays
    % Time is arbitrary
    
    % Equation 2.18
    pDelay      = rx.pDelay; 
    atten       = rx.atten;
    h           = zeros(1,length(t)); 
    for ii= 1:length(t)
        temp  = sum(atten(round(pDelay,10) ==  round(t(ii),10)));
        h(ii) = temp;
    end
    channelImp = resample(h, linspace(pDelay(1), pDelay(end), length(h))); % Might want to make this not start at zero
    dt = (length(channelImp)/(pDelay(end) - pDelay(1)))^(-1); 
    Nz = zeros(1, N); 
    channelImp = [0.*(0:dt:(pDelay(1)-dt)), channelImp, 0.*((pDelay - dt):dt:N*dt)]; 
    channelImpt = linspace(0, N*dt, length(channelImp)); 
            % Necessary for basic fourier theory and abstracting to
            % function of time...? 
    
end

function channelFreq = channelFrequencyResponse(freq, time)
    % I think both inputs should be an array

end


% function out = dirac2(tau, atten)
%     % Outputs an anonymous function that becomes 1 at tau
%     % This seems really dumb
%     out = @(t) []; 
%     for i = 1:length(tau)
%         temp = @(t) atten.*(round(tau, 4) == round(t, 4));
%         out = @(t) out + temp; 
%     end
% end

function h = dtbIR(recieverStruct)
    W = 1E7; fc = 1E8; % I
    M = 0:200; l = linspace(-4*pi, 4*pi); 

    h = zeros(1, length(l)); 
    for it = 1:length(recieverStruct.atten)
        ai  = recieverStruct.atten(it); 
        Ti  = recieverStruct.pDelay(it);
        aib = ai.*exp(-sqrt(-1).*2*pi*fc*Ti); 
        
            toBeSummed = aib.*sinc(l - Ti.*W); 
            if length(recieverStruct.atten) < 20
                yyaxis left; hold on; plot(abs(toBeSummed));  
                yyaxis right; hold off; plot(abs(h)); pause(0.1)
            end
            h = h + toBeSummed; % yyaxis right; plot(abs(h)); 
        % plot(abs(toBeSummed)); drawnow; pause(0.1); 
    end
    title('Effect of each scatterer'); ylabel('Sum'); 
    yyaxis left; ylabel('Individual'); 
    figure;
end
function rSt = delaySpread(rSt)
    W = 1E7; fc = 1E8; % I should not have to do this... 
    pDelays = rSt.pDelay;
    rSt.Td = range(pDelays); 
    rSt.Wc = 1/2/rSt.Td; 
    if      rSt.Wc > W/10; rSt.fading = 'flat';
    elseif  rSt.Wc < 10*W;  end
end
function rSt = dopplerSpread(rSt, Rx)
    % This seems like its much more complex than my assumptions...
    % Stationary transmitter and reciever really simplify this....
    
    % Bistatic sorta situation
    %   range-rate Tx-Ref --> 0; 
    %   Only need Ref-rec. Roving reciever
    W = 1E7; fc = 1E9; C = 3E8; 
    refPos = rSt.pos;
    for i = 1:length(Rx.pos) % Time iterater
        recPosX = Rx.pos(i, 1); recPosY = Rx.pos(i, 2);
        vecIn  = [recPosX- refPos(:, 1), recPosY - refPos(:, 2)]; 
        % Assume only the direct path from the reflector matters
        % I may be totally wrong with this algebra
        for ii = 1:length(vecIn)
            vecIn(ii, :)  = vecIn(ii,:)./norm(vecIn(ii,:));
            rangeR(ii,:)  = dot(vecIn(ii,:), Rx.V);
        end
        delayRate    = rangeR/C; 
        rSt.di(i, :) = fc.*      delayRate; 
        rSt.Ds(i)    = fc *range(delayRate); 
        rSt.Tc(i)   = 1/4/rSt.Ds(i);

    end
    


end
    

    
    
    
    
