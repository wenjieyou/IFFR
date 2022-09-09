% matlab实现并行批处理功能, 调用psiblast计算pssm矩阵
% wenjie (2017.08.26)

clear, clc
% % 启动 Matlab 并行计算环境
% CoreNum=12;     %设定机器CPU核心数量，CoreNum=12
% if matlabpool('size')<=0    %判断并行环境是否已启动
%     matlabpool('open','local',CoreNum);     %若未启动则启动并行环境
% else
%     disp('Already initialized');    %说明并行环境已经启动。
% end

d = dir('PDB1075_seq_FASTA');
seq_num = size(d,1)-2;

tic,
str1 = 'psiblast -in_msa PDB1075_seq_FASTA\seq_FASTA';
% str1 = 'psiblast -query PDB1075_seq_FASTA\seq_FASTA';         % 可计算DNA序列
str2 = '.fasta -db db_PDB1075 -out  PDB1075_pssm\PSSM';
str3 = '.txt -evalue 0.001 -num_iterations 3 -comp_based_stats 0 -out_ascii_pssm PDB1075_pssmMatrix\PSSMmatrix';
str4 = '.txt';

parfor i=1: seq_num
    cmd_str = [str1, num2str(i), str2, num2str(i), str3, num2str(i), str4];
    status = system(cmd_str);
    if status
        i
    end
end

toc

