%582 PCA Part 2 
clear all;
close all;
clc;

%% Camera 1
camera_1_data = load('cam1_2.mat');

camera_1_data.vidFrames1_2 = double(camera_1_data.vidFrames1_2);
[l w d numFrames1_2] = size(camera_1_data.vidFrames1_2);

for k = 1:numFrames1_2
    mov1_2(k).cdata = uint8(camera_1_data.vidFrames1_2(:,:,:,k));
    mov1_2(k).colormap = [];
end

for j = 1:numFrames1_2
   X_1_2 = rgb2gray(frame2im(mov1_2(j)));
   X1(:,:,j) = double(X_1_2);
   %  imshow(X_1_2); drawnow;
end


%% Camera 2
camera_2_data = load('cam2_2.mat');
camera_2_data.vidFrames2_2 = double(camera_2_data.vidFrames2_2);
[l w d numFrames2_2] = size(camera_2_data.vidFrames2_2);

for k = 1:numFrames2_2
    mov2_2(k).cdata = uint8(camera_2_data.vidFrames2_2(:,:,:,k));
    mov2_2(k).colormap = [];
end
for j = 1:numFrames2_2
     X_2_2 = rgb2gray(frame2im(mov2_2(j)));
   X2(:,:,j) = double(X_2_2);
   %  imshow(X_2_2); drawnow;
end


%% Camera 3
camera_3_data = load('cam3_2.mat');
camera_3_data.vidFrames3_2 = double(camera_3_data.vidFrames3_2);
[l w d numFrames3_2] = size(camera_3_data.vidFrames3_2);

for k = 1:numFrames3_2
    mov3_2(k).cdata = uint8(camera_3_data.vidFrames3_2(:,:,:,k));
    mov3_2(k).colormap = [];
end
for j = 1:numFrames3_2
    X_3_2  = rgb2gray(frame2im(mov3_2(j)));
   X3(:,:,j) = double(X_3_2);
  %  imshow(X_3_2); drawnow;
end

%% Data Extraction

XCrop = X1(250:480,:,:);
for n = 1:numFrames1_2
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
%      %   scatter(Acol1(n),Arow1(n));
%        scatter(col,row);
%          axis([0 640 0 480]);
%        pause(.05)
% %       hold on;
end
for n = 1:numFrames1_2
    if n ==1
    [row col] = find(X2(:,:,n)>=245); %find white spots
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
for n = 1:numFrames1_2
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
%        pause(.05)
%        grid on;
%       hold on;
end
% %NaN removal
Arow_2(isnan(Arow_2))=0;
Acol_2(isnan(Acol_2))=0;
%linar interpolation
for i = 1:length(Arow_2)
    if Arow_2(i)==0
        Arow_2(i) = (Arow_2(i-1) + Arow_2(i+1))/2;
        Acol_2(i) = (Acol_2(i-1) + Acol_2(i+1))/2;
    end
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
legend('Mode 1', 'Mode 2');
subplot(2,1,2);
plot(Y(:,1));
hold on
plot(Y(:,2));
legend('Mode 1', 'Mode 2');
%Modal Energy
figure(3)
semilogy(s/sum(s),'o','MarkerFaceColor',[0 0.447 0.741]);
grid on
title('Modal Energy');