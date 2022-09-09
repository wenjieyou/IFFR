% matlabʵ�ֲ�����������, ����psiblast����pssm����
% wenjie (2017.08.26)

clear, clc
% % ���� Matlab ���м��㻷��
% CoreNum=12;     %�趨����CPU����������CoreNum=12
% if matlabpool('size')<=0    %�жϲ��л����Ƿ�������
%     matlabpool('open','local',CoreNum);     %��δ�������������л���
% else
%     disp('Already initialized');    %˵�����л����Ѿ�������
% end

d = dir('PDB1075_seq_FASTA');
seq_num = size(d,1)-2;

tic,
str1 = 'psiblast -in_msa PDB1075_seq_FASTA\seq_FASTA';
% str1 = 'psiblast -query PDB1075_seq_FASTA\seq_FASTA';         % �ɼ���DNA����
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

