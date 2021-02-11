% Part 3 Comfortably Numb Guitar
clear all; clc
[y2, Fs2] = audioread('Floyd.m4a');
time2 = length(y2)/Fs2; % record time in seconds
L=10; n=length(y2)/time2*L;
t2=linspace(0,L,n+1); t=t2(1:n);
k=(1/L)*[0:n/2-1 -n/2:-1]; ks=fftshift(k);
% Bandpass filter
y_filter = bandpass(y2,[200,700],Fs2);
ygt_spec2=[];
tslide=0:0.03:t(end);
for j=1:length(tslide)
g=0.1*exp(-450*(t-tslide(j)).^2)'; % Gabor window
yg2=g.*y_filter(1:n); %GNR in Gabor Window time domain
ygt2=fft(yg2)'; %fft of gabor window GNR
ygt_spec2=[ygt_spec2; fftshift(abs(ygt2))];
end


figure(3)
pcolor(tslide,ks,log(abs(ygt_spec2).'+1));
shading interp
set(gca,'Fontsize',[14]);
colormap(hot)
axis([0 L 200 750])
xlabel('Time (t)','Fontsize',14), ylabel('frequency (Hz)','Fontsize',14);
title('Comfortably Numb Guitar','FontSize',18);
hold on, yline(293.66,'w','D4','Fontsize',16), hold on, yline(329.63,'w','E4','Fontsize',16),hold on
yline(369.99,'w','F#4','Fontsize',16), hold on, yline(392,'w','G4','Fontsize',16), hold on
yline(415.3,'w','G#4','Fontsize',16),hold on 
yline(440,'w','A4','Fontsize',16),hold on, yline(493.88,'w','B4','Fontsize',16), hold on
yline(554.37,'w','C#5','Fontsize',16),hold on, yline(587.33,'w','D5','Fontsize',16); 


