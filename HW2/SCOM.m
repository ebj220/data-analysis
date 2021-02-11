%%Part 1 Sweet Child'O Mine
[y1, Fs1] = audioread('GNR.m4a');
time1 = length(y1)/Fs1; % record time in seconds

L=5; %domain size (could be 10 seconds, hours, etc)
n=length(y1)/time1*L;
t2=linspace(0,L,n+1); 
t=t2(1:n);
k=(1/L)*[0:n/2-1 -n/2:-1]; %frequency Hz
ks=fftshift(k);

%% Gabor Filter
ygt_spec1=[]; 
tslide=0:0.05:t(end);
for j=1:length(tslide)
g=.05*exp(-5000*(t-tslide(j)).^2); % Gabor window size
yg1=g.*transpose(y1(1:n)); 
ygt1=fft(yg1); 
ygt_spec1=[ygt_spec1; fftshift(abs(ygt1))]; 
end

figure(1)
pcolor(tslide,ks,log(abs(ygt_spec1).'+1));
shading interp
set(gca,'Fontsize',[20]);
xlabel('Time (s)','FontSize',14), ylabel('Frequency (Hz)','FontSize',14),title('Sweet Child Of Mine Guitar','FontSize',18)
axis([0 L 0 1000]);
colormap(hot);
hold on, yline(277.18,'w','C#4','Fontsize',16), hold on, yline(311.13,'w','D#4','Fontsize',16),hold on
yline(369.99,'w','F#4','Fontsize',16), hold on, yline(415.3,'w','G#4','Fontsize',16), hold on, yline(554.37,'w','C#5','Fontsize',16);
hold on, yline(698.46,'w','F5','Fontsize',16), hold on, yline(739.99,'w','F#5','Fontsize',16);