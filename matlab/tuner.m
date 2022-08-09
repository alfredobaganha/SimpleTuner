clear all;
close all;
% load audio file
[y,Fs] = audioread('..\sample\guitar_sample.mp3');

% convert stereo to mono if that is the case
s = size(y);
if(s(2)~= 1)
    ym = mean(y,2); % couldn't thing of a better option
    clear s, y;
else
    ym = y;
    clear s, y;
end



% constant factors
scale = 1/100;
nHarm = 5;
nFft =16384;
filter_flag = false;

if(filter_flag)
    a = [1];
    b = Hlp.Numerator';
    ym = filter(b,a,ym);
end

% output vectors
acumm_LA = zeros(nFft,1);
acumm_Z = zeros(nFft,1);
fnv = zeros(nFft,nHarm);

for n=0:100

    if ~exist('att', 'var')
        att = 1;
        note_idx = [];
    else
        att = att+1;
    end
    
    startPoint = floor((length(ym)-scale*Fs)*rand());
    startPint = 875274+floor(n*441/2);

    yt = ym(startPoint:startPoint+scale*Fs-1);

    for i = 1:nHarm
        fnv(:,i) = fft(yt(1:floor(scale*Fs/i)).*hann(floor(scale*Fs/i)),nFft);
        acumm_Z = acumm_Z + fnv(:,i);
        acumm_LA = acumm_LA + 10*log10(abs(fnv(:,i)));
    end

    figure(1);
    hold on;
    grid on;
    la_acumm_Z = 10*log10(abs(acumm_Z));
    plot(la_acumm_Z(1:end/2));
    plot((acumm_LA(1:end/2)));
    %legend(['Acumm_Z ';'Acumm_LA']);

    [mxz, mxz_i] = max(la_acumm_Z(1:end/2));
    [mxl, mxl_i] = max(acumm_LA(1:end/2));
    plot([mxz_i, mxl_i],[mxz, mxl],'rx')

    note_idx(att)=mxz_i;

end

disp(note_idx');