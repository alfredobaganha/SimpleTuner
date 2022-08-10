clear all;
close all;
% load audio file
[y,Fs] = audioread('..\sample\guitar_sample.mp3'); % sample file
%[y,Fs] = audioread('d:\castle.mp3');

% convert stereo to mono if that is the case
s = size(y);
if(s(2)~= 1)
    ym = mean(y,2); % couldn't think of a better solution
    clear s, y;
else
    ym = y;
    clear s, y;
end


% Not explicit anywhere, but we assume standart A4 = 440Hz and equal
% temperament

% constant factors
scale = 1;                  % scale factor

nHarm = 5;                  % number of harmonic components. 5 seemed to 
                            % work consistently, 3 although simpler tended
                            % to have octave resolution problems

nFft =65536;                % fft size @SR=44.1kHz (cd quality) frequency 
                            % resolution is about 0.67Hz. this is far more
                            % than we actually need but help the low
                            % frequency spectrum, particularly C#1 (9 string
                            % guitar) at the cost of expensive computing
                            % and higher frequncy error
                       
filter_flag = false;        % enable/disable optional FIR LFP

if(filter_flag)
    a = [1];
    load('filter.mat');     % LPF with 2.5kHz cutoff and 4kHz stop 40dB rejection
    ym = filter(b,a,ym);
end

figure(1);
title('Complex combination');

for n=0:10 % lets check about 100ms of the input sample

    % output vectors
    %acumm_LA = zeros(nFft,1);
    mul_Z = ones(nFft,1);
    fnv = zeros(nFft,nHarm);
    
    if ~exist('att', 'var')
        att = 1;
        note_idxz = [];
        note_idxl = [];
    else
        att = att+1;
    end
    
    if(n==0)
        startPoint = floor((length(ym)-scale*Fs)*rand());   % randomizing starting point
    else
        startPoint = startPoint + n*scale*Fs/100;           %advance 10ms
    end
    %startPoint = 1+floor(n*nFft/100);
    %startPoint = 875274+floor(n*441);

    yt = ym(startPoint:startPoint+scale*Fs-1);        % extract scale*1s worth of samples

    for i = 1:nHarm
        if(i==1)
		inp = yt.*hann(length(yt));     % This is the fundamental pitch.
                                        % hanning windowing to mitigate 
                                        % high frequency noise
                                        % TODO: Consider more suitable
                                        % window function
        else
            inp = resample(inp,1,i);    % downsample with an polyphase anti-alias filter
                                        % to remove higher order harmonics
        end
        fnv(:,i) = fft(inp,nFft);       % the trick here is to zero pad the fft.
                                        % storing is not necessarily
                                        % required
        
		mul_Z = mul_Z.*fnv(:,i);        % Harmonic Product Spectrum: The idea
                                        % is to coherently multiply the
                                        % multiple filtered inputs. the
                                        % fundamental harmonic will combine
                                        % coherently while higher orde
                                        % harmonics not so much, thus given
                                        % a distintive peak at the
                                        % funtamental fequency pitch
    end

   
    la_mul_Z = 10*log10(abs(mul_Z));
    
    figure(1);
    %subplot(2,1,1);
    hold on; grid on; xlim([50,1960]);
    plot(la_mul_Z(1:end/2));
    [mxz, mxz_i] = max(la_mul_Z(1:end/2));
    plot(mxz_i,mxz,'rx');
%     subplot(2,1,2);
%     hold on; grid on;%xlim([14,491]);
%     plot(unwrap(angle(mul_Z(1:end/2))));

    note_idxz(att)=mxz_i - 1; % -1 to correct for matlab's not being zero indexed

end

disp([note_idxz']); % check spreadsheet for the index/note correspondence +- 1 error is äcceptable