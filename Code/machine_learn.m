  %HW4 Classifying Digits
clc; clear all; close all;

%Load in data
[images_train, labels_train] = mnist_parse('train-images-idx3-ubyte', 'train-labels-idx1-ubyte');
[images_test, labels_test] = mnist_parse('t10k-images-idx3-ubyte', 't10k-labels-idx1-ubyte');

for i = 1:6000
    X_train(:,i) = reshape(images_train(:,:,i),[784 1]);
end
for i = 1:1000
    X_test(:,i) = reshape(images_test(:,:,i+2000),[784 1]);
end
X_train = double(X_train);
X_test = double(X_test);
labels_train = labels_train(1:6000);
labels_test = labels_test(2001:3000);

%% Sorting
[labels_train_sort, Indextrain] = sort(labels_train,'ascend');
X_train_sort=X_train(:,Indextrain);
[labels_test_sort, Indextest] = sort(labels_test,'ascend');
X_test_sort=X_test(:,Indextest);
%% SVD 
[U, S, V] = svd(X_train_sort,'econ');
[Us, Ss, Vs] = svd(X_test_sort,'econ');
Xr_test = Ss(1:250,1:250)*Vs(:,1:250)';
Xr_train = S(1:250,1:250)*V(:,1:250)';
figure(1)
subplot(2,2,1)
 plot(diag(S)/sum(diag(S)),'ro','Linewidth',[2])
 ylabel('Modal Energy');
 xlabel('Number of Modes');
%% Modal Truncation 
for j = 1:1
    Xm1 = U(:,1:j)*S(1:j,1:j)*V(:,1:j)';
end
for j = 1:10
    Xm10 = U(:,1:j)*S(1:j,1:j)*V(:,1:j)';
end
for j = 1:50
    Xm50 = U(:,1:j)*S(1:j,1:j)*V(:,1:j)';
end
for j = 1:250
    Xm250 = U(:,1:j)*S(1:j,1:j)*V(:,1:j)';
end
for j = 1:350
    Xm350 = U(:,1:j)*S(1:j,1:j)*V(:,1:j)';
end
figure(2)
subplot(2,3,1);
imshow(reshape(X_train_sort(:,600),[28 28]));
title('Orignial Image');
subplot(2,3,2);
imshow(reshape(Xm1(:,600),[28 28]));
title('Single Mode Reconstruction');
subplot(2,3,3);
imshow(reshape(Xm10(:,600),[28 28]));
title('10 Mode Reconstruction')
subplot(2,3,4);
imshow(reshape(Xm50(:,600),[28 28]));
title('50 Mode Reconstruction');
subplot(2,3,5);
imshow(reshape(Xm250(:,600),[28 28]));
title('250 Mode Reconstruction');
subplot(2,3,6);
imshow(reshape(Xm350(:,600),[28 28]));
title('350 Mode Reconstruction');
%% Count Data
c1 = 0; c2 = 0; c3 = 0; c4 = 0; c5 = 0; c6 = 0; 
c7 = 0; c8 = 0; c9 = 0; c0 = 0;
sort_data = zeros
for i = 1:length(labels_train)
    if labels_train(i) == 1
        c1 = c1+1;
    end
    if labels_train(i) == 2
        c2 = c2+1;
    end
    if labels_train(i) == 3
        c3 = c3+1;
    end
    if labels_train(i) == 4
        c4 = c4+1;
    end
    if labels_train(i) == 5
        c5 = c5+1;
    end
    if labels_train(i) == 6
        c6 = c6+1;
    end
    if labels_train(i) == 7
        c7 = c7+1;
    end
    if labels_train(i) == 8
        c8 = c8+1;
    end
    if labels_train(i) == 9
        c9 = c9+1;
    end
    if labels_train(i) == 0
        c0 = c0+1;
    end
end

count1 = 0; count2 = 0; count3 = 0; count4 = 0; count5 = 0; count6 = 0; 
count7 = 0; count8 = 0; count9 = 0; count0 = 0;
sort_data = zeros
for i = 1:length(labels_test)
    if labels_test(i) == 1
        count1 = count1+1;
    end
    if labels_test(i) == 2
        count2 = count2+1;
    end
    if labels_test(i) == 3
        count3 = count3+1;
    end
    if labels_test(i) == 4
        count4 = count4+1;
    end
    if labels_test(i) == 5
        count5 = count5+1;
    end
    if labels_test(i) == 6
        count6 = count6+1;
    end
    if labels_test(i) == 7
        count7 = count7+1;
    end
    if labels_test(i) == 8
        count8 = count8+1;
    end
    if labels_test(i) == 9
        count9 = count9+1;
    end
    if labels_test(i) == 0
        count0 = count0+1;
    end
end
zero = [1:c0];
one = [c0+1:c0+c1];
two = [c0+c1+1:c0+c1+c2];
three = [c0+c1+c2+1:c0+c1+c2+c3];
four = [c0+c1+c2+c3+1:c0+c1+c2+c3+c4];
five = [c0+c1+c2+c3+c4+1:c0+c1+c2+c3+c4+c5];
six = [c0+c1+c2+c3+c4+c5+1:c0+c1+c2+c3+c4+c5+c6];
seven = [c0+c1+c2+c3+c4+c5+c6+1:c0+c1+c2+c3+c4+c5+c6+c7];
eight = [c0+c1+c2+c3+c4+c5+c6+c7+1:c0+c1+c2+c3+c4+c5+c6+c7+c8];
nine = [c0+c1+c2+c3+c4+c5+c6+c7+c8+1:c0+c1+c2+c3+c4+c5+c6+c7+c8+c9];
%% Grouping
Zero = Xr_test(:,1:85);
One = Xr_test(:,86:211);
Two = Xr_test(:,212:327);
Three = Xr_test(:,328:434);
Four = Xr_test(:,435:544);
Five = Xr_test(:,545:631);
Six = Xr_test(:,632:718);
Seven = Xr_test(:,719:817);
Eight = Xr_test(:,818:906);
Nine = Xr_test(:,907:end);

Zero_t = Xr_train(:,1:c0);
One_t = Xr_train(:,one);
Two_t = Xr_train(:,two);
Three_t = Xr_train(:,three);
Four_t = Xr_train(:,four);
Five_t = Xr_train(:,five);
Six_t = Xr_train(:,six);
Seven_t = Xr_train(:,seven);
Eight_t = Xr_train(:,eight);
Nine_t = Xr_train(:,nine);
for i = length(labels_train_sort)
    if labels_train_sort(i) == 0
        Zero_t = Xr_train(:,i);
    end
    if labels_train_sort(i) == 1
        One_t = Xr_train(:,i);
    end
end

%% Problem 4

 figure(3)
 %subplot(2,2,1);
 scatter3(V(zero,2),V(zero,3),V(zero,5));
 %plot zeros
 hold on
 scatter3(V(one,2),V(one,3),V(one,5),'filled');
 %plot ones
 hold on
scatter3(V(two,2),V(two,3),V(two,5),'filled');
%plot twos
 hold on
scatter3(V(three,2),V(three,3),V(three,5),'filled');
%plot threes
 hold on
 scatter3(V(four,2),V(four,3),V(four,5),'filled');
%plot fours
 hold on
 scatter3(V(five,2),V(five,3),V(five,5),'filled');
%plot fives
 hold on
 scatter3(V(six,2),V(six,3),V(six,5),'go','filled');
% %plot sixes
 hold on
scatter3(V(seven,2),V(seven,3),V(seven,5),'ro','filled');
%plot sevens
 hold on
scatter3(V(eight,2),V(eight,3),V(eight,5),'mo','filled');
%plot eights
 hold on
 scatter3(V(nine,2),V(nine,3),V(nine,5),'yo','filled');
%plot nines
legend('0','1','2','3','4','5','6','7','8','9');
title('PCA Projected Data Clusters');

[idx0,C0] = kmeans([V(zero,2) V(zero,3) V(zero,5)],1);
[idx1,C1] = kmeans([V(one,2) V(one,3) V(one,5)],1);
[idx2,C2] = kmeans([V(two,2) V(two,3) V(two,5)],1);
[idx3,C3] = kmeans([V(three,2) V(three,3) V(three,5)],1);
[idx4,C4] = kmeans([V(four,2) V(four,3) V(four,5)],1);
[idx5,C5] = kmeans([V(five,2) V(five,3) V(five,5)],1);
[idx6,C6] = kmeans([V(six,2) V(six,3) V(six,5)],1);
[idx7,C7] = kmeans([V(seven,2) V(seven,3) V(seven,5)],1);
[idx8,C8] = kmeans([V(eight,2) V(eight,3) V(eight,5)],1);
[idx9,C9] = kmeans([V(nine,2) V(nine,3) V(nine,5)],1);

figure(4);
%subplot(2,2,2);
plot3(C0(:,1),C0(:,2),C0(:,3),'x','Markersize',15,'LineWidth',3);
hold on
plot3(C1(:,1),C1(:,2),C1(:,3),'x','Markersize',15,'LineWidth',3);
hold on
plot3(C2(:,1),C2(:,2),C2(:,3),'x','Markersize',15,'LineWidth',3);
hold on
plot3(C3(:,1),C3(:,2),C3(:,3),'x','Markersize',15,'LineWidth',3);
hold on
plot3(C4(:,1),C4(:,2),C4(:,3),'x','Markersize',15,'LineWidth',3);
hold on
plot3(C5(:,1),C5(:,2),C5(:,3),'x','Markersize',15,'LineWidth',3);
hold on
plot3(C6(:,1),C6(:,2),C6(:,3),'gx','Markersize',15,'LineWidth',3);
hold on
plot3(C7(:,1),C7(:,2),C7(:,3),'rx','Markersize',15,'LineWidth',3);
hold on
plot3(C8(:,1),C8(:,2),C8(:,3),'mx','Markersize',15,'LineWidth',3);
hold on
plot3(C9(:,1),C9(:,2),C9(:,3),'yx','Markersize',15,'LineWidth',3);
grid on;
legend('0','1','2','3','4','5','6','7','8','9');
title('K-Means Centroids');
%% Linear Discriminant Analysis
%initialize pairs and triplets of data
Test_Pair = [Four Nine];
Train_Pair = [Four_t Nine_t];
label_Train_Pair = [4*ones(c4,1); 9*ones(c9,1)];
label_Test_Pair = [4*ones(count4,1); 9*ones(count9,1)];

%LDA classification line
pre2 = classify(transpose(Test_Pair),transpose(Train_Pair),label_Train_Pair);
counter = 0;
for i =1:length(pre2)
    if pre2(i) == label_Test_Pair(i);
        counter = counter+1;
    end
end
fprintf('LDA Percent Error is %4.2f %', ((length(pre2)-counter)/length(pre2))*100);
figure(5)
scatter(1:length(pre2),pre2,'ro','filled');
ylabel('Digit Number');
axis([1 length(pre2) 0 9]);
% % 
%% SVM classifier with training data, labels and test set
  Mdl = fitcecoc(transpose(X_train_sort),labels_train_sort);
  pre3 = predict(Mdl,transpose(X_test_sort));
  
  counter1 = 0;
for i =1:length(pre3)
    if pre3(i) == labels_test_sort(i);
        counter1 = counter1+1;
    end
end
fprintf('SVM Percent Error is %4.2f %', ((length(pre3)-counter1)/length(pre3))*100);
figure(6)
scatter(1:length(pre3),pre3,'bo','filled');
ylabel('Digit Number');
axis([1 length(pre3) 0 9]);

% %% Classification Tree

 tree=fitctree(transpose(X_train_sort),labels_train_sort);
% view(tree,'mode','graph');
 pre4 = predict(tree,transpose(X_test_sort));

 %data plotting and error calcuation
 counter4 = 0;
for i =1:length(pre4)
    if pre4(i) == labels_test_sort(i);
        counter4 = counter4+1;
    end
end
fprintf('Classification Tree Percent Error is %4.2f %', ((length(pre4)-counter4)/length(pre4))*100);
figure(7)
scatter(1:length(pre4),pre4,'go','filled');
ylabel('Digit Number');
axis([1 length(pre4) 0 9]);