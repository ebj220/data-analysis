%DMD Image Analysis
clc; clear all; close all;
%% Load Videos
ski_data = VideoReader('ski_drop_low.mp4');

vidFrames = read(ski_data);
numFrames = get(ski_data,'numberOfFrames');
for k = 1 : numFrames
mov(k).cdata = vidFrames(:,:,:,k);
mov(k).colormap = [];
end

for j=1:numFrames
Xg=rgb2gray(frame2im(mov(j)));
Xski(:,:,j) = double(Xg);
%imshow(Xg); drawnow
end
%% DMD analysis F1 Video
%defining spaciotemporal domain
t = linspace(0 ,6,numFrames);
dt = t(2)-t(1);
x = linspace(1,51800);

%initializing data matrix
for j = 1:numFrames
Xdat(:,j) = reshape(Xski(:,:,j),[518400 1]);
end
%creation of DMD data matrices
Xa = Xdat(:,1:end-1);
Xb = Xdat(:,2:end);
%SVD of Xa
[U,S,V] = svd(Xa,'econ');

figure(1)
plot(diag(S)/sum(diag(S)),'ro','Linewidth',[2]) 
xlabel('SVD Modes')
ylabel('Modal Energy Content')
r = 30;
%rank truncation of SVD matricies
Ur = U(:,1:r);
Sr = S(1:r,1:r);
Vr = V(:,1:r); 

Atilde = Ur'*Xb*Vr*inv(Sr); 
[W,D] = eig(Atilde);
Phi = Xb*Vr*inv(Sr)*W; %DMD modes (eigenvectors)

lambda = diag(D); %eigenvalues
omega = log(lambda)/dt; %frequencies

%plot omega to determine background modes
figure(2)
plot(omega,'.','LineWidth',7);
hold on
xline(0);
yline(0);
xlabel('Real')
ylabel('Complex')
%% Background Subtracktion
%find indicies of low value frequencies
Index_b = find(abs(omega)<.01);
Ingex_f = setdiff(1:r, Index_b);

omega_b = omega(Index_b);
Phi_b = Phi(:,Index_b);
B = Phi_b\Xdat(:,1);

for t2 = 1:length(t)
    time_dynamics(:,t2) = B.*exp(omega_b.* t(t2));
end
X_dmd = Phi_b*time_dynamics;
%Background image reconstruction
X_for = Xdat-abs(X_dmd);
%Forground image reconstruction

figure(3)
for j = 1:3
    frame2 = reshape(X_for(:,100 +(j-1)*40),[540 960]);
    frame2 = uint8(140-real(frame2));
    subplot(2,3,j)
    imshow(frame2);
    title(['Frame ', num2str(100+(j-1)*40)])
    axis('equal','tight')
end
for j = 1:3
    frame1 = reshape(X_dmd(:,100 +(j-1)*40),[540 960]);
    frame1 = uint8(real(frame1));
    subplot(2,3,3+j)
    imshow(frame1);
    title(['Frame ', num2str(100+(j-1)*40)])
    axis('equal','tight')
end
