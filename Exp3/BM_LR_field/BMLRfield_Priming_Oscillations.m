%% write by mengjin;

%% 8conds: {1=LFF; 2=RFF; 3=LHH; 4=RHH;
%%          5=LHF; 6=RHF; 7=LFH; 8=RFH;}
clear; clc;
filename = fullfile(pwd,"BMLRfield_raw_detrnd_data.mat");
load(filename);

testLR = RT_detrnd;

subjnum = size(testLR, 1);
condnum = size(testLR, 2);
ISInum  = size(testLR, 3);

testLR_diff = zeros(subjnum, condnum/2, ISInum);
for i = 1: condnum/2
    testLR_diff(:,i,:) = testLR(:,i+4,:)-testLR(:,i,:);
end


%% set FFT parameters
fs = 1/0.02;  % 50Hz
fft_length = 2^7;
hann_w = hann(fft_length);
zeropadding = zeros(fft_length - ISInum, 1);

for subj_idx = 1 : subjnum

    for cond_idx = 1: condnum

        temp = squeeze(testLR(subj_idx, cond_idx, :));

        if ~iscolumn(temp)
            temp = temp';
        end

        % ?跺～??浠ュ?归?? FFT ?垮害
        s = [temp; zeropadding];

        s = s.* hann_w;

        [Y(subj_idx,cond_idx, :), f] = plotfft(s, fft_length, fs);
    end

    for condiff_idx = 1:condnum/2

        tempdiff = squeeze(testLR_diff(subj_idx, condiff_idx, :));

        s_diff = [tempdiff; zeropadding];

        s_diff = s_diff.* hann_w;

        [Y_diff(subj_idx, condiff_idx, :), f] = plotfft(s_diff, fft_length, fs);
    end
end

amp_Y = squeeze(mean(abs(Y), 1));% mean across subjs, return 8*128 matrix
amp_Y_diff = squeeze(mean(abs(Y_diff), 1));

%% 计算相位差

angle_all = angle(Y);
for i = 1:condnum/2
    angle_diff(:, i, :) = squeeze(angle(Y(:,i+4,:))-angle(Y(:,i,:)));
end

% 计算相位的余弦和正弦平方和通常与相位差的概率密度函数（PDF）有关
% cos^2(θ) + sin^2(θ)
cv_angle_diff = squeeze(mean(cos(angle_diff)).^2 + mean(sin(angle_diff)).^2);


% 如果是计算相位一致性，一个更常见的度量是相位锁定值（Phase Locking Value, PLV），
% 它通常用于评估两个信号在特定频率上的相位同步性。
% PLV的计算方法如下：
complex_phases = exp(1i * angle_all);
for cond_idx = 1:condnum/2
    for subj_idx = 1:subjnum
        angle_cond1 = complex_phases(subj_idx, cond_idx,:);
        angle_cond2 = complex_phases(subj_idx, cond_idx+4,:);
        plv(subj_idx, :) = squeeze(angle_cond1 .* conj(angle_cond2));
    end
    Plvs(cond_idx, :) = abs(mean(plv, 1)); % mean across subjs
end

%% bootstraping -- 置换检验

permutation_size = 1000;

% 初始化变量
testLR_permu = zeros(subjnum, condnum, ISInum);

amp_Y_permu = zeros(permutation_size, condnum, fft_length);
amp_Y_diffpermu = zeros(permutation_size, condnum/2, fft_length);

for k = 1: permutation_size

    Y_permu = zeros(subjnum, condnum, fft_length);
    Y_diff_permu = zeros(subjnum, condnum/2, fft_length);

    for subj_idx = 1:subjnum
        %% FFT for all conds permu_data

        for cond_idx = 1:condnum

            rand_idx = randperm(ISInum); % shuffing the ISI dimension

            s_permu = squeeze(testLR(subj_idx,cond_idx,rand_idx));
            testLR_permu(subj_idx,cond_idx,:) = s_permu;
            
            % 确保s_permu是列向量            
            s_permu = columnize(s_permu);
            
            % 零填充和应用Hann窗
            s_permu = [s_permu; zeropadding];
            s_permu = s_permu .* hann_w;

            % 计算FFT
            [Y_permu(subj_idx,cond_idx,:),f] = plotfft(s_permu,fft_length,fs);
        end
        
        %% FFT for conds difference (IC-C) permu_data
        for idx = 1:4

            s_diff_permu = squeeze(testLR_permu(subj_idx, idx+4, :) ...
                    - testLR_permu(subj_idx, idx, :));

            s_diff_permu = columnize(s_diff_permu);

            % 零填充和应用Hann窗
            s_diff_permu = [s_diff_permu; zeropadding];
            s_diff_permu = s_diff_permu .* hann_w;
            
            % 计算FFT
            [Y_diff_permu(subj_idx,idx,:),f] = plotfft(s_diff_permu,fft_length,fs);
        
            % 计算FFT之后IC与C之间的相位差
            angle_diff_permu(subj_idx,idx,:) = squeeze(angle(Y_permu(subj_idx,idx+4,:)) ...
                         - angle(Y_permu(subj_idx, idx, :)));

        end
    end

    % 计算每次排列后的幅度
    amp_Y_permu(k,:,:) = squeeze(mean(abs(Y_permu), 1)); % mean across subjects
    amp_Y_diffpermu(k,:,:) = squeeze(mean(abs(Y_diff_permu), 1));

    cv_angle_diff_permu(k, :, :) = mean(cos(angle_diff_permu),1).^2 ...
        + mean(sin(angle_diff_permu),1).^2; 
    %  return= permutation size * 4 * fftlength

end


F_range = find((f >= 0) & (f <= 12)); % frequency index
alpha = 0.005;
freq = 1;

T_amp = zeros(condnum, fft_length);
T_amp_diff = zeros(condnum, fft_length);
for n = 1:condnum
    data = amp_Y(n, F_range)';
    permu_data = squeeze(amp_Y_permu(:, n, F_range));
    [sig, test_threshold_H, threshold_H] = permutest(data,permu_data,alpha,freq);

    T_amp(n, :) = test_threshold_H.*ones(size(f)); 
%     T_amp(n, F_range) = threshold_H;

end

data=[]; permu_data=[];
for n = 1:4
    data = amp_Y_diff(n, F_range)';
    permu_data = squeeze(amp_Y_diffpermu(:, n, F_range));
    [sig_diff, test_threshold_H_diff, threshold_H_diff] = permutest(data,permu_data,alpha,freq);

    T_amp_diff(n, :) = test_threshold_H_diff.*ones(size(f)); 
%     T_amp_diff(n, F_range) = threshold_H_diff;
end



%% %% figure(1); frequency and power of 8 conditions %%%%
titles = {'target face in the LVF', ...
          'target face in the RVF',...
          'target house in the LVF',...
          'target house in the RVF'};
lengeds = {'Left Face-Face (LFF)',   'Left House-Face (LHF)'; ...
           'Right Face-Face (RFF)',  'Right House-Face (RHF)'; ...
           'Left House-House (LHH)',  'Left Face-House (LFH)'; ...
           'Right House-House (RHH)', 'Right Face-House (RFH)'};

figure;
for i = 1:4  
    subplot(2,2,i)
    hold on; plot(f,amp_Y(i, :),'color','r','linewidth',3);
    hold on; plot(f,amp_Y(i+4, :),'color','k','linewidth',3); 
    hold on; plot(f, T_amp(i, :),'--r','linewidth',1); % thershold
    hold on; plot(f, T_amp(i+4, :),'--k','linewidth',1); 
    set(gca,'xlim',[0 12]); 
    set(gca,'xtick', [0:3:12]);
    set(gca,'ylim',[0 0.3]); 
    set(gca,'ytick', [0:0.3:0.3]);
    xlabel('Frequency'); 
    ylabel('Amplitude');
    title(sprintf('FFT for %s', titles{i})); %, 'Position', [0.5, 1.05, 0]
    legend(lengeds{i, 1}, lengeds{i, 2}, 'Box', 'off', 'Color', 'none');
end
hold off;
set(gcf, 'Color', 'none');
filename = 'results_BMLRfield all conditions oscillations.png';
% saveas(gcf, filename);

%% %% figure(2) ; frequency and power of 4(C-IC) condition 

legends = {'LFF-LHF', 'RFF-RHF', 'LHH-LFH', 'RHH-RFH'};
figure;
for i = 1:4  
    subplot(2,2,i); hold on; 
    plot(f, amp_Y_diff(i, :),'color','k', 'LineStyle', '-', 'linewidth', 3); 
    plot(f, amp_Y_diff(i, :),'color','r', 'LineStyle', '--','linewidth', 3); 
    plot(f, T_amp_diff(i, :),'color','r','LineStyle', '--','linewidth',2); % thershold
    plot(f, T_amp_diff(i, :),'color','k','LineStyle', '--','LineWidth',2); 
    set(gca,'xlim',[0 12]); 
    set(gca,'xtick', [0:3:12]);
    set(gca,'ylim',[0 0.3]); 
    set(gca,'ytick', [0:0.3:0.3]);
    xlabel('Frequency'); 
    ylabel('Amplitude');
    title(sprintf('FFT for %s', titles{i})); 
    legend(legends{i}, 'Box', 'off', 'Color', 'none');
end
hold off;
set(gcf, 'Color', 'none');
filename = 'results_BMLRfield IC-C oscillation.png';
% saveas(gcf, filename);







