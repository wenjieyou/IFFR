% ��ǰ�����psiblast�õ���pssm���󣬴�����mat��ʽ
% ��ȡָ��Ŀ¼�����е�PSSMmatrix�ļ�
% wenjie (2021.09.17)

clear, clc
tic
filename1='PSSMmatrix';
filename2='.txt';
dirname='PhosP_pssmMatrix_nr2020';
d = dir(dirname);
fn = size(d,1)-2; % ��Ŀ¼���ļ���(��ȥ��ǰĿ¼����Ŀ¼2��)
pssm = cell(fn,1);
for k=1:fn
    k
    % header �������ж����
    [header,sequence]=fastaread(['PhosP_FASTA\seq_FASTA', num2str(k), '.fasta']);
    seqlen=size(sequence,2);   %��ȡ���г���
    filename=[filename1, num2str(k), filename2];
    tmp1=importdata([dirname, '\', filename]);
    % ��ȡ��Ŀ¼��PSSMmatrix?.txt�е���ֵ
    dat=tmp1.data;
    blastPSSM = dat(1:seqlen,1:20);       % pssm����seqlen��*20��
%     blastPSSM = 1./(1+exp(-blastPSSM));   % sigmoid��һ��
    pssm{k} = blastPSSM;
end

label=[ones(1132,1); 0*ones(638,1)];
toc
save PhosP_nr2020_1132_638_pssmMatrix pssm label