%582 PCA 
clear all;
clc;
close all;
%% Camera 1
camera_1_data = load('cam1_3.mat');

camera_1_data.vidFrames1_3 = double(camera_1_data.vidFrames1_3);
[l w d numFrames1_3] = size(camera_1_data.vidFrames1_3);

for k = 1:numFrames1_3
    mov1_3(k).cdata = uint8(camera_1_data.vidFrames1_3(:,:,:,k));
    mov1_3(k).colormap = [];
end

for j = 1:numFrames1_3
   X_1_3 = rgb2gray(frame2im(mov1_3(j)));
   X1(:,:,j) = double(X_1_3);
   %  imshow(X_1_3); drawnow;
end
%% Camera 2
camera_2_data = load('cam2_3.mat');
camera_2_data.vidFrames2_3 = double(camera_2_data.vidFrames2_3);
[l w d numFrames2_3] = size(camera_2_data.vidFrames2_3);

for k = 1:numFrames2_3
    mov2_3(k).cdata = uint8(camera_2_data.vidFrames2_3(:,:,:,k));
    mov2_3(k).colormap = [];
end
for j = 1:numFrames2_3
     X_2_3 = rgb2gray(frame2im(mov2_3(j)));
   X2(:,:,j) = double(X_2_3);
    % imshow(X_2_3); drawnow;
end

%% Camera 3
camera_3_data = load('cam3_3.mat');
camera_3_data.vidFrames3_3 = double(camera_3_data.vidFrames3_3);
[l w d numFrames3_3] = size(camera_3_data.vidFrames3_3);

for k = 1:numFrames3_3
    mov3_3(k).cdata = uint8(camera_3_data.vidFrames3_3(:,:,:,k));
    mov3_3(k).colormap = [];
end
for j = 1:numFrames3_3
    X_3_3 = rgb2gray(frame2im(mov3_3(j)));
   X3(:,:,j) = double(X_3_3);
   % imshow(X_3_3); drawnow;
end
%% Location Data Collection

%Camera 1
XCrop = X1(250:480,:,:);
for n = 1:numFrames3_3
    if n ==1
    [row col] = find(XCrop(:,:,n)>=240); %find white spots
    for i = 1:length(row)
        row(i) = row(i) + 250;
    end
    Arow_1(n) = sum(row)/length(row);
    Acol_1(n) = sum(col)/length(col);
    else
        [row col]=find(X1(:,:,n)>=240);
        i = 1;
        while i <= length(row)
            if row(i)<(Arow_1(n-1)-30)|row(i)>(Arow_1(n-1)+30)
                row(i)=[];
                col(i)=[];
            else
                i = i+1;
            end
        end
        i = 1;
         while i <= length(row)
            if col(i)<(Acol_1(n-1)-30)|col(i)>(Acol_1(n-1)+30)
                row(i)=[];
                col(i)=[];
            else
                i = i+1;
            end
        end
        Arow_1(n) = sum(row)/length(row);
        Acol_1(n) = sum(col)/length(col);
    end
%         
% %      %   scatter(Acol1(n),Arow1(n));
%         scatter(col,row);
%         axis([0 640 0 480]);
%         pause(.05)
% %       hold on;
end
%Camera 2
for n = 1:numFrames3_3
    if n ==1
    [row col] = find(X2(:,:,n)>=255); %find white spots
    Arow_2(n) = sum(row)/length(row);
    Acol_2(n) = sum(col)/length(col);
    else
        [row col]=find(X2(:,:,n)>=247);
        i = 1;
        while i <= length(row)
            if row(i)<(Arow_2(n-1)-30)|row(i)>(Arow_2(n-1)+30)
                row(i)=[];
                col(i)=[];
            else
                i = i+1;
            end
        end
        i = 1;
         while i <= length(row)
            if col(i)<(Acol_2(n-1)-20)|col(i)>(Acol_2(n-1)+20)
                row(i)=[];
                col(i)=[];
            else
                i = i+1;
            end
        end
        Arow_2(n) = sum(row)/length(row);
        Acol_2(n) = sum(col)/length(col);
    end
        
%      %   scatter(Acol1(n),Arow1(n));
%        scatter(col,row);
%        axis([0 640 0 480]);
%        pause(.05)
%        grid on;
%       hold on;
end
% Camera3

for n = 1:numFrames3_3
    if n ==1
    [row col] = find(X3(:,:,n)>=255); %find white spots
    Arow_3(n) = sum(row)/length(row);
    Acol_3(n) = sum(col)/length(col);
    else
        [row col]=find(X3(:,:,n)>=235);
        i = 1;
        while i <= length(row)
            if row(i)<(Arow_3(n-1)-45)|row(i)>(Arow_3(n-1)+45)
                row(i)=[];
                col(i)=[];
            else
                i = i+1;
            end
        end
        i = 1;
         while i <= length(row)
            if col(i)<(Acol_3(n-1)-45)|col(i)>(Acol_3(n-1)+45)
                row(i)=[];
                col(i)=[];
            else
                i = i+1;
            end
        end
        Arow_3(n) = sum(row)/length(row);
        Acol_3(n) = sum(col)/length(col);
    end
        
%      %   scatter(Acol1(n),Arow1(n));
%        scatter(col,row);
%        axis([0 640 0 480]);
%        pause(0.05);
%       grid on;
%       hold on;
end
 %% PCA Analysis
X = [Acol_1;Arow_1;Acol_2;Arow_2;Acol_3;Arow_3]; %data matrix
[m,n]=size(X); %size of X
Xave = mean(X,2);
Xm=X-repmat(Xave,1,n); %subtract mean of data
% scatter3(Xm(1,:),Xm(2,:),Xm(3,:));

[u,s,v]=svd(Xm'/sqrt(n-1),'econ'); % perform the SVD
Y=u.*X'; % produce the principal components projection
figure(2)
subplot(2,1,1);
plot(u(:,1));
hold on
plot(u(:,2));
hold on
plot(u(:,3));
legend('Mode 1', 'Mode 2','Mode 3');
subplot(2,1,2);
plot(Y(:,1));
hold on
plot(Y(:,2));
hold on
plot(Y(:,3));
legend('Mode 1', 'Mode 2','Mode 3');
%Modal Energy
figure(3)
semilogy(s/sum(s),'o','MarkerFaceColor',[0 0.447 0.741]);
grid on
title('Modal Energy');