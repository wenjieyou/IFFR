% matlab实现并行批处理功能, 调用psiblast计算pssm矩阵
% wenjie (2022.10.04)

clear, clc

d = dir('PDB1075_seq_FASTA');
seq_num = size(d,1)-2;

tic,
str1 = 'psiblast -query PDB1075_seq_FASTA\seq_FASTA';         % 可计算DNA序列
str2 = '.fasta -db uniref50 -out  PDB1075_pssm\PSSM';
str3 = '.txt -evalue 0.001 -num_iterations 3 -comp_based_stats 0 -out_ascii_pssm PDB1075_pssmMatrix\PSSMmatrix';
str4 = '.txt';

parfor i=1: seq_num
    fn_pssmMatrix = ['PDB1075_pssmMatrix\PSSMmatrix', num2str(i), '.txt'];
    if exist(fn_pssmMatrix, 'file')==0      % 文件不存在
        cmd_str = [str1, num2str(i), str2, num2str(i), str3, num2str(i), str4];
        status = system(cmd_str);
        if status
            i
        end
    end
end

toc

