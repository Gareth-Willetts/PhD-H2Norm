% p = [5 4 3 2 1];
% q = [6 5 4 3 2];
% 
% tic
% pZeros = length(q) + 1;
% qZeros = length(p) + 1;
% 
% pp = [p zeros(1, pZeros)];
% length(pp)
% qq = [q zeros(1, qZeros)];
% length(qq)
% ifft(fft(pp) .* fft(qq))
% toc
% 
% tic
% conv(p,q)
% toc
% 
% %%
% timeFFT_means = [];
% timeConv_means = [];
% for i = 1000:1020
%     timeFFT = [];
%     timeConv = [];
%     for j = 1:10000
%         p = 1:i;
%         q = 2:i+1;
%         tic
%         convfft(p,q);
%         % pZeros = length(q) + 1;
%         % qZeros = length(p) + 1;
%         % 
%         % pp = [p zeros(1, pZeros)];
%         % %length(pp)
%         % qq = [q zeros(1, qZeros)];
%         % %length(qq)
%         % ifft(fft(pp) .* fft(qq));
%         timeFFT(i-4) = toc;
% 
%         tic
%         conv(p,q);
%         timeConv(i-4) = toc;
%     end
%     timeFFT_means(i-999) = mean(timeFFT);
%     timeConv_means(i-999) = mean(timeConv);
% end
% 
% figure;
% plot(timeFFT_means); hold on; plot(timeConv_means);
% legend("FFT", "conv")

%%
close all; clear all; clc;
timeFFT_means = [];
timeConv_means = [];
for i = 5:2000
    timeFFT = [];
    timeConv = [];
    for j = 1:10000
        p = 1:i;
        q = 2:i+1;
        tic
        convnfft(p,q);
        % pZeros = length(q) + 1;
        % qZeros = length(p) + 1;
        % 
        % pp = [p zeros(1, pZeros)];
        % %length(pp)
        % qq = [q zeros(1, qZeros)];
        % %length(qq)
        % ifft(fft(pp) .* fft(qq));
        timeFFT(i-4) = toc;
    
        tic
        conv(p,q);
        timeConv(i-4) = toc;
    end
    timeFFT_means(i-4) = mean(timeFFT);
    timeConv_means(i-4) = mean(timeConv);
end

%%
x = 5:2000;
figure;
plot(x(1:16),timeFFT_means(1:16)); hold on; plot(x(1:16),timeConv_means(1:16));
xlabel("Polynomial Degree"); ylabel("Time (s)"); title("Polynomial Multiplication Comparison");
set(gca, 'FontSize', 22)
legend("FFT", "conv")

%%
x = 5:2000;
figure;
plot(x(996:end),timeFFT_means(996:end)); hold on; plot(x(996:end),timeConv_means(996:end));
xlabel("Polynomial Degree"); ylabel("Time (s)"); title("Polynomial Multiplication Comparison");
set(gca, 'FontSize', 22)
legend("FFT", "conv")