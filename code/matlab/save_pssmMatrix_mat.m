% 对前面计算psiblast得到的pssm矩阵，存贮成mat格式
% 读取指定目录下所有的PSSMmatrix文件
% wenjie (2021.09.17)

clear, clc
tic
filename1='PSSMmatrix';
filename2='.txt';
dirname='PhosP_pssmMatrix_nr2020';
d = dir(dirname);
fn = size(d,1)-2; % 子目录下文件数(减去当前目录、父目录2个)
pssm = cell(fn,1);
for k=1:fn
    k
    % header 可用于判断类别
    [header,sequence]=fastaread(['PhosP_FASTA\seq_FASTA', num2str(k), '.fasta']);
    seqlen=size(sequence,2);   %读取序列长度
    filename=[filename1, num2str(k), filename2];
    tmp1=importdata([dirname, '\', filename]);
    % 读取子目录下PSSMmatrix?.txt中的数值
    dat=tmp1.data;
    blastPSSM = dat(1:seqlen,1:20);       % pssm矩阵seqlen行*20列
%     blastPSSM = 1./(1+exp(-blastPSSM));   % sigmoid归一化
    pssm{k} = blastPSSM;
end

label=[ones(1132,1); 0*ones(638,1)];
toc
save PhosP_nr2020_1132_638_pssmMatrix pssm label