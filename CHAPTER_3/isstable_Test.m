%%
close all; clear all; clc;
% %% Execute for 5, 15, 25, 35, 45
% % Create vectors for storing results
% syms s
% counter = 1;
% time_inbuilt_means = zeros(1,5);
% time_poly_means = zeros(1,5);
% time_poly_numerical_means = zeros(1,5);
% numIters = 100000;
% 
% for j = 5:10:45   
%     time_inbuilt = zeros(1,numIters);
%     time_poly = zeros(1,numIters);
%     time_poly_numerical = zeros(1,numIters);
%     symb = (s+1)^(j);
%     cd = sym2poly(symb);
%     cn = ones(1,length(cd)-1);
%     for i = 1:numIters
%         tic
%         inbuiltResult = isstable(tf(cn,cd));
%         time_inbuilt(i) = toc;
%         tic
%         polyResult = PolynomialStabilityTest(cd,cn);
%         time_poly(i) = toc;
%         tic
%         polyNumericalResult = PolynomialNumericalStabilityTest(cd,cn);
%         time_poly_numerical(i) = toc;
%     end
%     time_inbuilt_means(counter) = mean(time_inbuilt(10:end));
%     time_poly_means(counter) = mean(time_poly(10:end));
%     time_poly_numerical_means(counter) = mean(time_poly_numerical(10:end));
%     counter = counter + 1;
% end

%% Execute for entire range (5-45)
% syms s
% counter = 1;
% time_inbuilt_means = zeros(1,41);
% time_poly_means = zeros(1,41);
% time_poly_numerical_means = zeros(1,41);
% numIters = 100000;
% 
% for j = 5:1:45   
%     time_inbuilt = zeros(1,numIters);
%     time_poly = zeros(1,numIters);
%     time_poly_numerical = zeros(1,numIters);
%     symb = (s+1)^(j);
%     cd = sym2poly(symb);
%     cn = ones(1,length(cd)-1);
%     for i = 1:numIters
%         tic
%         inbuiltResult = isstable(tf(cn,cd));
%         time_inbuilt(i) = toc;
%         tic
%         polyResult = PolynomialStabilityTest(cd,cn);
%         time_poly(i) = toc;
%         tic
%         polyNumericalResult = PolynomialNumericalStabilityTest(cd,cn);
%         time_poly_numerical(i) = toc;
%     end
%     time_inbuilt_means(counter) = mean(time_inbuilt(10:end));
%     time_poly_means(counter) = mean(time_poly(10:end));
%     time_poly_numerical_means(counter) = mean(time_poly_numerical(10:end));
%     counter = counter + 1;
% end
% 
% %%
% figure;
% plot(5:45, time_inbuilt_means, 'r')
% hold on;
% plot(5:45, time_poly_means, 'b')
% plot(5:45, time_poly_numerical_means, 'g')
% legend("isstable", "Symbolic Polynomial Stability Test", "Numerical Polynomial Stability Test", 'Location', 'northwest')
% xlabel("n"); ylabel("Time (s)");
% title("Comparison of stability test times")
% set(gca, "FontSize", 20)

%% Larger range (only numerical version)
% syms s
% counter = 1;
% time_inbuilt_means = zeros(1,41);
% time_poly_numerical_means = zeros(1,41);
% numIters = 100000;
% 
% for j = 5:1:85  
%     time_inbuilt = zeros(1,numIters);
%     time_poly_numerical = zeros(1,numIters);
%     symb = (s+1)^(j);
%     cd = sym2poly(symb);
%     cn = ones(1,length(cd)-1);
%     for i = 1:numIters
%         tic
%         inbuiltResult = isstable(tf(cn,cd));
%         time_inbuilt(i) = toc;
%         tic
%         polyNumericalResult = PolynomialNumericalStabilityTest(cd,cn);
%         time_poly_numerical(i) = toc;
%     end
%     time_inbuilt_means(counter) = mean(time_inbuilt(10:end));
%     time_poly_numerical_means(counter) = mean(time_poly_numerical(10:end));
%     counter = counter + 1;
% end
% 
% %%
% figure;
% plot(5:85, time_inbuilt_means, 'r')
% hold on;
% plot(5:85, time_poly_numerical_means, 'g')
% legend("isstable", "Numerical Polynomial Stability Test", 'Location', 'northwest')
% xlabel("n"); ylabel("Time (s)");
% title("Comparison of stability test times")
% set(gca, "FontSize", 20)

%%
% load("isstable_results_UpTo45.mat")
% figure;
% plot(5:45, time_inbuilt_means, 'r')
% hold on;
% plot(5:45, time_poly_means, 'b')
% plot(5:45, time_poly_numerical_means, 'g')
% legend("isstable", "Symbolic Polynomial Stability Test", "Numerical Polynomial Stability Test", 'Location', 'northwest')
% xlabel("n"); ylabel("Time (s)");
% title("Comparison of stability test times")
% set(gca, "FontSize", 20)
% 
% time_inbuilt_means(1:10:41)

%%
load("isstable_results_UpTo85.mat")
figure;
plot(5:85, time_inbuilt_means, 'r')
hold on;
plot(5:85, time_poly_numerical_means, 'g')
legend("isstable", "Numerical Polynomial Stability Test", 'Location', 'northwest')
xlabel("n"); ylabel("Time (s)");
title("Comparison of stability test times")
set(gca, "FontSize", 20)