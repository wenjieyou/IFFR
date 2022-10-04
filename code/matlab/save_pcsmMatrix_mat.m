% ���㲢�������������е�����ѧ���ԣ���������������PSSM��ʽ��
% ���������ں�PSSM���ж���Ϣ�ںϼ��㡣
% wenjie (2022.10.04)

clear, clc
tic,
load LiuBin6_KBS

pcdat = liubin6;

aa_str=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'];
pcsm=cell(1075,1);
for i=1:1075
    fidin=fopen(['PDB1075_BPs_525_550_pssmMatrix\PSSMmatrix',num2str(i), '.txt'],'r');
    k=1;
    id=zeros(1,1);
    while ~feof(fidin)
        tline=fgetl(fidin);
        
        % ��ֵ�����е��ص㣺�к�(��ֵ)��AA(��ĸ)���÷�(��ֵ)��...��
        tmp=isletter(tline);        % ���������ص��tmp��ֻ��1��1������
        if sum(tmp)==1
            res=tline(tmp);
            id(k)=findstr(res,aa_str);
            k = k+1;
        end
    end
    fclose(fidin);
    
    % header �������ж����
    [header,sequence]=fastaread(['PDB1075_seq_FASTA\seq_FASTA', num2str(i), '.fasta']);
    seqlen=size(sequence,2);   %��ȡ���г���
    
    temp=pcdat(id,:);
    pcsm{i}=temp(1:seqlen,:);
end

save PDB1075_LiuBin6_pcsm6Matrix pcsm
toc