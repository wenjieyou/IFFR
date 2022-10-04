% 计算并存贮蛋白质序列的物理化学属性，并存贮成类似于PSSM格式，
% 该数据用于和PSSM进行多信息融合计算。
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
        
        % 数值存贮行的特点：行号(数值)，AA(字母)，得分(数值)，...，
        tmp=isletter(tline);        % 符合上面特点的tmp是只含1个1的向量
        if sum(tmp)==1
            res=tline(tmp);
            id(k)=findstr(res,aa_str);
            k = k+1;
        end
    end
    fclose(fidin);
    
    % header 可用于判断类别
    [header,sequence]=fastaread(['PDB1075_seq_FASTA\seq_FASTA', num2str(i), '.fasta']);
    seqlen=size(sequence,2);   %读取序列长度
    
    temp=pcdat(id,:);
    pcsm{i}=temp(1:seqlen,:);
end

save PDB1075_LiuBin6_pcsm6Matrix pcsm
toc