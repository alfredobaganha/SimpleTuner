clear all;
close all;
clc;
% Totally unecessary

note_names = ['A_0____';'A#/Bb_0';'B_0____';'C_1____';'C#/Db_1';'D_1____';...
              'D#/Eb_1';'E_1____';'F_1____';'F#/Gb_1';'G_1____';'G#/Ab_1';...
              'A_1____';'A#/Bb_1';'B_1____';'C_2____';'C#/Db_2';'D_2____';...
              'D#/Eb_2';'E_2____';'F_2____';'F#/Gb_2';'G_2____';'G#/Ab_2';...
              'A_2____';'A#/Bb_2';'B_2____';'C_3____';'C#/Db_3';'D_3____';...
              'D#/Eb_3';'E_3____';'F_3____';'F#/Gb_3';'G_3____';'G#/Ab_3';...
              'A_3____';'A#/Bb_3';'B_3____';'C_4____';'C#/Db_4';'D_4____';...
              'D#/Eb_4';'E_4____';'F_4____';'F#/Gb_4';'G_4____';'G#/Ab_4';...
              'A_4____';'A#/Bb_4';'B_4____';'C_5____';'C#/Db_5';'D_5____';...
              'D#/Eb_5';'E_5____';'F_5____';'F#/Gb_5';'G_5____';'G#/Ab_5';...
              'A_5____';'A#/Bb_5';'B_5____';'C_6____';'C#/Db_6';'D_6____';...
              'D#/Eb_6';'E_6____';'F_6____';'F#/Gb_6';'G_6____';'G#/Ab_6';...
              'A_6____';'A#/Bb_6';'B_6____';'C_7____';'C#/Db_7';'D_7____';...
              'D#/Eb_7';'E_7____';'F_7____';'F#/Gb_7';'G_7____';'G#/Ab_7';...
              'A_7____';'A#/Bb_7';'B_7____';'C_8____'];

% load audio file
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
baseV = [2; 3; 5; 6; 7];
expV = (fix(log(Fs)./log(baseV))+1);
res = (Fs./(baseV.^expV));
opt = mean(res);
[~,idx] = min(abs(res -opt));
nFft = baseV(idx)^expV(idx);
clear baseV expV res opt idx;


% Not explicit anywhere, but we assume equal temperament
std_pitch = 440;            % Hz - A4 reference pitch
n_ratio = power(2,1/12);
scale = (std_pitch/16).*(n_ratio).^(0:87)'; % Piano scale from A0 to C8 (88 keys)
extremes = [scale(1)/n_ratio ; scale(88)*n_ratio];
% constant factors
nHarm = 5;                  % number of harmonic components. 5 seemed to 
                            % work consistently, 3 although simpler tended
                            % to have octave resolution problems


                       
filter_flag = false;        % enable/disable optional FIR LFP
plot_flag = false;

if(filter_flag)             % Assume Fs = 44100
    a = [1];
    load('filter.mat');     % LPF with 2.5kHz cutoff and 4kHz stop 40dB rejection
    ym = filter(b,a,ym);
end

windowSize = max(1,fix(Fs/20)); % about 50ms
step = max(1,fix(Fs/100));      % about 10ms
n_max = floor((length(ym) - 2*windowSize)/step);

%n_max = 10;

if(plot_flag)
    figure(1);
    title('Harmonic Product Spectrum');
end

note_idx = zeros(n_max,1);
t_adv =zeros(n_max,1);

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
    [mxz, mxz_i] = max(la_mul_Z(1:fix(nFft/2)));
    
    if(plot_flag)
        figure(1);
        hold on; grid on; %xlim([50,1960]);
        plot(la_mul_Z(1:fix(nFft/2)));
        plot(mxz_i,mxz,'rx');
    end

    note_idx(att)=mxz_i - 1; % -1 to correct for matlab's not being zero indexed
    t_adv(att) = startPoint-1;
    
end

ID_note_idx = note_idx;
clear note_idxz;
ID_note_f = ID_note_idx*Fs/nFft;

% limit scope to +- half step from the considered scale
t = (ID_note_f>=extremes(1))&(ID_note_f<=extremes(2));
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

