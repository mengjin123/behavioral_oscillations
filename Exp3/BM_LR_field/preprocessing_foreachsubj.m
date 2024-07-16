function testLR = preprocessing_foreachsubj(ISI_list, n_conds)
% å¯¹æ?°æ??????å§?RTè¿?è¡?ç¬?ä¸?æ¬¡æ?æ´?ï¼?
% ä¸»è????????¤é??è¯?trialï¼?ä»¥å??RT??2.5ä¸?????å·?ä¹?å¤???ï¼?
% testLR size = subj * cond(8 conds) * ISI *repeats (32 trials(per conds))
% 
%%

data = load(fullfile(pwd, "BMLRfield_allsubjects_data.mat"));
subjnum = size(data.all_subjects_data, 1);
n_isi   = length(ISI_list);

allsubjs_RTzscore = zeros(subjnum, n_conds, n_isi, 32);
figure; hold on;
for subj_index = 1:subjnum
     
    a = squeeze(data.all_subjects_data(subj_index,1));
    a = cell2mat(a);
    n_trials = size(a, 1);
    n_repeatTrials = n_trials/(n_conds * n_isi);
    
    %% ?¾å??Response(=0)å¯¹å???è¡?ç´¢å?ï¼?å¹¶å?è¿?äº?è¡?å¯¹å???RT=NANï¼?
    accuracy = a(:, 9);
    RT_cond = a(:,8);
      
    RT_cond(accuracy==0) = NaN; % find(accuracy==0); index1=a(:, 9)== 0; 
    RT_cond(RT_cond==0)  = NaN;
    
    Upper_threshold = nanmean(RT_cond) + 2 * nanstd(RT_cond);
    Lower_threshold = nanmean(RT_cond) - 2 * nanstd(RT_cond);
    RT_cond(RT_cond >= Upper_threshold) = NaN;
    RT_cond(RT_cond <= Lower_threshold) = NaN;
    
    accuracy(isnan(RT_cond))=0;
    Correct_response = length(find(accuracy(:,1)==1))./n_trials;
    
    mean_val = nanmean(RT_cond); % mean(RT_cond(~isnan(RT_cond)));
    meidan_val = nanmedian(RT_cond);
    RT_stds = nanstd(RT_cond);
    fprintf("\n subj%d: \n accuracy = %.02f, RT_std= %.02f, \n " + ...
            "meanRT = %.02f, medianRT = %.02f, \n " + ...
            "maxRT = %.02f, minRT = %.02f \n ", ...
            subj_index, ...
            Correct_response, RT_stds,...
            mean_val, meidan_val, ...
            max(RT_cond), min(RT_cond));
    
    
    %% zscore
    RT_zscore = (RT_cond - nanmean(RT_cond))/nanstd(RT_cond);

    %% plot each subj RT distribution
    subplot(3, 4, subj_index);
    h = histogram(RT_cond, round(RT_stds*100));
    set(h, 'FaceAlpha', 0.5, 'EdgeColor', 'black');
    titles = sprintf(['S%d: std = %.02f, \n ' ...
                    'meanRT = %.02f,'],...
                subj_index, RT_stds, mean_val);
    title(titles);

    %% 8conds: {1=LFF; 2=RFF; 3=LHH; 4=RHH; 
    %%          5=LHF; 6=RHF; 7=LFH; 8=RFH;}
    b(:,1) = (a(:,4)==1)&(a(:,3)==1)&(a(:,6)==1)&(a(:,5)==1); % LFLF
    b(:,2) = (a(:,4)==2)&(a(:,3)==1)&(a(:,6)==2)&(a(:,5)==1); % RFRF
    b(:,3) = (a(:,4)==1)&(a(:,3)==2)&(a(:,6)==1)&(a(:,5)==2); % LHLH
    b(:,4) = (a(:,4)==2)&(a(:,3)==2)&(a(:,6)==2)&(a(:,5)==2); % RHRH

    b(:,5) = (a(:,4)==1)&(a(:,3)==2)&(a(:,6)==1)&(a(:,5)==1); % LHLF
    b(:,6) = (a(:,4)==2)&(a(:,3)==2)&(a(:,6)==2)&(a(:,5)==1); % RHRF
    b(:,7) = (a(:,4)==1)&(a(:,3)==1)&(a(:,6)==1)&(a(:,5)==2); % LFLH
    b(:,8) = (a(:,4)==2)&(a(:,3)==1)&(a(:,6)==2)&(a(:,5)==2); % RFRH

    
    FG1 = nan(n_conds, n_isi, n_repeatTrials);
    for i = 1:n_conds    
        for k = 1:length(ISI_list)
            index = zeros(n_repeatTrials,1);
            index(:) = find((abs(a(:,2)-ISI_list(k))<=1e-3)& b(:,i)== 1);
            FG1(i, k, :) = RT_zscore(index);
        end
    end
    allsubjs_RTzscore(subj_index, :, :, :) = FG1; 
    % subj * cond(8 conds) * ISI *repeats (32 trials(per conds))
end
hold off;

xlabel('Reaction Time (RT)');
ylabel('Frequency');

filename = 'BMLRfield_RTdistribution_allsubjs.png';
saveas(gcf, filename);

testLR = allsubjs_RTzscore;
end





