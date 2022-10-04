% ��ǰ�����psiblast�õ���pssm���󣬴�����mat��ʽ
% ��ȡָ��Ŀ¼�����е�PSSMmatrix�ļ�
% wenjie (2022.10.04)

clear, clc
tic
filename1='PSSMmatrix';
filename2='.txt';
dirname='PDB1075_pssmMatrix_nr2020';
d = dir(dirname);
fn = size(d,1)-2; % ��Ŀ¼���ļ���(��ȥ��ǰĿ¼����Ŀ¼2��)
pssm = cell(fn,1);
for k=1:fn
    k
    % header �������ж����
    [header,sequence]=fastaread(['PDB1075_FASTA\seq_FASTA', num2str(k), '.fasta']);
    seqlen=size(sequence,2);   %��ȡ���г���
    filename=[filename1, num2str(k), filename2];
    tmp1=importdata([dirname, '\', filename]);
    % ��ȡ��Ŀ¼��PSSMmatrix?.txt�е���ֵ
    dat=tmp1.data;
    blastPSSM = dat(1:seqlen,1:20);       % pssm����seqlen��*20��
%     blastPSSM = 1./(1+exp(-blastPSSM));   % sigmoid��һ��
    pssm{k} = blastPSSM;
end

label=[ones(525,1); 0*ones(550,1)];
toc
save PDB1075_pssmMatrix pssm label