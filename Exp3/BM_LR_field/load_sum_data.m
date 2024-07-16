%% write by Mengjin Li, at 2024/6/4;

clc; 
clear;
all_subjects_data = cell(12, 1);

for subj_num = 2:13
    % 构造被试文件夹的名称
    subj_folder = sprintf('subj%d', subj_num);

    % 构造被试文件夹的完整路径
    subj_path = fullfile(subj_folder);

    % 初始化一个cell数组，用于存储当前被试的所有run数据
    subject_runs_data = cell(8, 1);

    for run_num = 1:8

        % 构造当前run的文件名
        file_name = sprintf('s%d_run%d_Matrix.mat', subj_num, run_num);

        % 构造完整的文件路径
        file_path = fullfile(subj_path, file_name);

        % 检查文件是否存在
        if exist(file_path, 'file') == 2 % 2表示检查文件是否存在

            % 加载.mat文件
            load(file_path,'Matrix', '-mat');

            % 将数据存储到相应的cell中
            subject_runs_data{run_num} = Matrix;
%             all_subjects_data{subj_num-1, run_num} = MatrixData;

        else
            % 如果文件不存在，记录一条警告，并用空数组填充cell
            warning('文件不存在：%s', file_path);
            subject_runs_data{run_num} = [];           
        end

        subject_data = cat(1, subject_runs_data{:});

    end

    %% for ConsciLR: 如果是第13个被试，只加载存在的run
%     if subj_num == 13 
%         empty_runs = arrayfun(@isempty,subject_runs_data);
%         if any(empty_runs)
%             subject_runs_data = subject_runs_data(~empty_runs);
%         end
%         
%     end
%     subject_data = cat(1, subject_runs_data{:});

    all_subjects_data{subj_num-1} = subject_data;
end

savefilename = 'BMLRfield_allsubjects_data.mat';
save(savefilename, 'all_subjects_data');







