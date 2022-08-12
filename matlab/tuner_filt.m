clear all;
close all;
clc;

% Totally unecessary ... But ...
note_names = ['A-0____';'A#/Bb-0';'B-0____';'C-1____';'C#/Db-1';'D-1____';...
              'D#/Eb-1';'E-1____';'F-1____';'F#/Gb_1';'G-1____';'G#/Ab-1';...
              'A-1____';'A#/Bb-1';'B-1____';'C-2____';'C#/Db-2';'D-2____';...
              'D#/Eb-2';'E-2____';'F-2____';'F#/Gb-2';'G-2____';'G#/Ab-2';...
              'A-2____';'A#/Bb-2';'B-2____';'C-3____';'C#/Db-3';'D-3____';...
              'D#/Eb-3';'E-3____';'F-3____';'F#/Gb-3';'G-3____';'G#/Ab-3';...
              'A-3____';'A#/Bb-3';'B-3____';'C-4____';'C#/Db-4';'D-4____';...
              'D#/Eb-4';'E-4____';'F-4____';'F#/Gb-4';'G-4____';'G#/Ab-4';...
              'A-4____';'A#/Bb-4';'B-4____';'C-5____';'C#/Db-5';'D-5____';...
              'D#/Eb-5';'E-5____';'F-5____';'F#/Gb-5';'G-5____';'G#/Ab-5';...
              'A-5____';'A#/Bb-5';'B-5____';'C-6____';'C#/Db_6';'D-6____';...
              'D#/Eb-6';'E-6____';'F-6____';'F#/Gb_6';'G-6____';'G#/Ab-6';...
              'A-6____';'A#/Bb-6';'B-6____';'C-7____';'C#/Db_7';'D_7____';...
              'D#/Eb-7';'E-7____';'F-7____';'F#/Gb-7';'G-7____';'G#/Ab-7';...
              'A-7____';'A#/Bb-7';'B-7____';'C-8____'];

% load sample audio file
[y,Fs] = audioread('..\sample\flute_a4_440.mp3'); % sample file
%[y,Fs] = audioread('..\sample\oboe_a4_440.mp3'); % sample file
%[y,Fs] = audioread('..\sample\piano_a4_440.mp3'); % sample file
%[y,Fs] = audioread('..\sample\violin_open_bow_pizz.mp3'); % sample file
%[y,Fs] = audioread('..\sample\violin_a.mp3'); % sample file
%[y,Fs] = audioread('..\sample\guitar_sample.mp3'); % sample file
[y,Fs] = audioread('d:\castle.mp3');

% convert stereo to mono if that is the case
s = size(y);
if(s(2)~= 1)
    ym = mean(y,2); % couldn't think of a better solution
    clear s;
else
    ym = y;
    clear s;
end

nFft =2^(fix(log2(Fs)+1));  % fft size @SR=44.1kHz (cd quality) frequency 
                            % resolution is about 0.67Hz. this is far more
                            % than we actually need but help the low
                            % frequency spectrum, particularly C#1 (9 string
                            % guitar) at the cost of expensive computing
                            % and higher frequncy error. nFft = 2^16

% Try to guess the optimum nFft size based on some simple "primes" that
% would achieve an optimum fft efficiency/frequency resolution.
% This is a tradeoff
baseV = [2; 3; 5; 6; 7; 10; 11];
expV = (fix(log(Fs)./log(baseV))+1);
res = (Fs./(baseV.^expV));
opt = mean(res);
[~,idx] = min(abs(res -opt));
nFft = baseV(idx)^expV(idx);

disp(['Sampling Rate Fs = ',num2str(Fs),'Hz']);
disp(['FFT: base = ',num2str(baseV(idx)),' expoent = ', num2str(expV(idx))]);
clear baseV expV res opt idx;

% Not explicit anywhere, but we assume equal temperament
std_pitch = 440;            % Hz - A4 reference pitch
disp(['Reference Pitch: A4 = ',num2str(std_pitch),'Hz']);

n_ratio = power(2,1/12);
scale = (std_pitch/16).*(n_ratio).^(0:87)'; % Piano scale from A0 to C8 (88 keys)
extremes = [scale(1)/n_ratio ; scale(88)*n_ratio];
disp('Temperament: Equal Temperament');

% constant factors
nHarm = 5;                  % number of harmonic components. 5 seemed to 
                            % work consistently, 3 although simpler tended
                            % to have octave resolution problems
                       
%filter_flag = false;        % enable/disable optional FIR LFP
plot_flag = false;

%[Deprecated]
% if(filter_flag)             % Assume Fs = 44100
%     a = [1];
%     load('filter.mat');     % LPF with 2.5kHz cutoff and 4kHz stop 40dB rejection
%     ym = filter(b,a,ym);
% end

windowSize = max(1,fix(Fs/10)); % about 100ms
step = max(1,fix(Fs/100));      % about 10ms
n_max = floor((length(ym) - 2*windowSize)/step);

%n_max = 10;

if(plot_flag)
    figure(1);
    title('Harmonic Product Spectrum');
end



note_idx = zeros(n_max,1);
t_adv =zeros(n_max,1);

load('fir_set');    % Loads the fir coefficients

for n=0:n_max % lets check about 100ms of the input sample

    if(mod(n,100)==0)
        str = sprintf('Working: %f%% complete',100*n/n_max);
        disp(str);
    end
    
    % output vectors
    mul_Z = ones(nFft,1);
    fnv = zeros(nFft,nHarm);
    
    if ~exist('att', 'var')
        att = 1;
    else
        att = att+1;
    end
    
    if(n==0)
        startPoint = 1;
    else
        startPoint = startPoint + step;                     % advance 10ms
    end
    

    yt = ym(startPoint:startPoint+windowSize-1);            % take 100ms of data
    yt = yt.*hann(length(yt));
    
    for i = 1:nHarm

        switch(i)
            case 1
                inp = yt;
                
            case 2
                inp = conv(yt, b2, 'same');

            case 3
                inp = conv(yt, b3, 'same');

            case 4
                inp = conv(yt, b4, 'same');

            case 5
                inp = conv(yt, b5, 'same');

            otherwise
                inp = yt;
        end

        fnv(:,i) = fft(inp, nFft);      % the trick here is to zero pad the fft.
                                        % storing is not necessarily
                                        % required
        
		mul_Z = mul_Z.*fnv(:, i);       % Harmonic Product Spectrum: The idea
                                        % is to coherently multiply the
                                        % multiple filtered inputs. the
                                        % fundamental harmonic will combine
                                        % coherently while higher orde
                                        % harmonics not so much, thus given
                                        % a distintive peak at the
                                        % funtamental fequency pitch
%         plot(10*log10(abs(mul_Z(1:round(nFft/2)))));
%         if(i==5)
%             pause(2);
%         end
    end

   
    la_mul_Z = 10*log10(abs(mul_Z));
    [mxz, mxz_i] = max(la_mul_Z(1:round(nFft/2)));
    
    if(plot_flag)
        figure(1);
        hold on; grid on; %xlim([50,1960]);
        plot(la_mul_Z(1:round(nFft/2)));
        plot(mxz_i,mxz,'rx');
    end

    note_idx(att)=mxz_i - 1; % -1 to correct for matlab's not being zero indexed
    t_adv(att) = startPoint-1;
    
end

ID_note_idx = note_idx;
clear note_idxz;
ID_note_f = ID_note_idx*Fs/nFft;

% limit scope to +- half step from the considered scale
t = (ID_note_f>extremes(1))&(ID_note_f<extremes(2));
ID_note_idx = ID_note_idx(t);
ID_note_f = ID_note_f(t);
t_adv = t_adv(t);
clear t;

% Identify the closest note in the scale
index = zeros(length(ID_note_f),1);
P_note_f = index;
note_err = index;                       % not compensating for the quantizaton error (yet)


for i=1:length(index)
    [~,index(i)]=min(abs(scale-ID_note_f(i)));
    P_note_f(i) = scale(index(i));
    note_err(i) = ID_note_f(i) - P_note_f(i);    
end

f = fopen('TunerLog.txt','w');

fprintf(f,'Fs = %dHz\nn_fft = %d\nReference Pitch A4 = %fHz\nEqual Temperament\n\n',Fs,nFft,std_pitch);
str = sprintf('T[ms]\tID_NT_f\tP_NT_f\tF_ERR\tNT_Name');
fprintf(f,'%s\n',str);
disp(str);

for i=1:length(index)
    str=sprintf('%fms\t%f\t%f\t%f\t%s',1000*t_adv(i)/Fs,ID_note_f(i),P_note_f(i),note_err(i),note_names(index(i),:));
    fprintf(f,'%s\n',str);
    disp(str);    
end

fclose(f);

