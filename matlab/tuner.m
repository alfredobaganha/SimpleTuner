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



 figure(1);
title('Complex combination');
hold on;
grid on;
xlim([14,491]);
figure(2)
title('Log-Lin combination');
hold on;
grid on;
xlim([14,491]);

for n=0:100

    % output vectors
    acumm_LA = zeros(nFft,1);
    acumm_Z = zeros(nFft,1);
    fnv = zeros(nFft,nHarm);
    
    if ~exist('att', 'var')
        att = 1;
        note_idxz = [];
        note_idxl = [];
    else
        att = att+1;
    end
    
    %startPoint = floor((length(ym)-scale*Fs)*rand());
    startPoint = 875274+floor(n*441/2);
    %startPoint = 875274+floor(n*441);

    yt = ym(startPoint:startPoint+scale*Fs-1);

    for i = 1:nHarm
        fnv(:,i) = (1/(nHarm -i +1))*fft(yt(1:floor(scale*Fs/i)).*hann(floor(scale*Fs/i)),nFft);
        acumm_Z = acumm_Z + fnv(:,i);
        acumm_LA = acumm_LA + 10*log10(abs(fnv(:,i)));
    end

   
    la_acumm_Z = 10*log10(abs(acumm_Z));
    figure(1);
    subplot(2,1,1);
    hold on; grid on; xlim([14,491]);
    plot(la_acumm_Z(1:end/2));
    [mxz, mxz_i] = max(la_acumm_Z(1:end/2));
    plot(mxz_i,mxz,'rx');
    subplot(2,1,2);
    hold on; grid on;xlim([14,491]);
    plot(unwrap(angle(acumm_Z(1:end/2))));
    
    
    figure(2);
    plot((acumm_LA(1:end/2)));
    [mxl, mxl_i] = max(acumm_LA(1:end/2));
    plot(mxl_i,mxl,'rx')
    %legend(['Acumm_Z ';'Acumm_LA']);

    note_idxz(att)=mxz_i;
    note_idxl(att)=mxl_i;

end

disp([note_idxz',note_idxl']);