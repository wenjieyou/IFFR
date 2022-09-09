% 袁超提供的蛋白质二级格式数据，计算PSSM
% 对python数据格式.npy的数据文件，转存为matlab数据格式.mat
%       from scipy import io
%       mat=np.load('CullPDB13114.npy')
%       io.savemat('seq.mat',{'seq':mat})
% 然后matlab把这些序列写入fasta格式
% wenjie (2081.11.24)

clear, clc
load seq.mat

for i=1:size(seq,1)
    tmp1 = seq(i,:,:);          % 三维的张量格式
    tmp2 = shiftdim(tmp1);      % 三维中抽取一个二维面
    tmp3 = tmp2(1,:);               % 所取出第一维，即原始序列
    tmp3(find(isspace(tmp3))) = [];     % 删除中间可能的空格
    len = length(tmp3);
    
    filename = ['seq_FASTA\seq_FASTA_', num2str(i)];
    fid = fopen(filename,'a');
    firstline = ['>', num2str(i)];          % 首先写入fasta格式：'>'
    fprintf(fid, '%s\n', firstline);
    
    int = fix(len/60);          %对序列长度取整数
    rem = mod(len, 60);     %对序列长度求余
    if int >1
        for j = 1:fix(len/60)
            tline = tmp3((j-1)*60+1:j*60);
            fprintf(fid, '%s\n', tline);
        end
    end
    if rem > 0
        fprintf(fid, '%s', tmp3(end+1-rem:end));
    end
    fclose(fid);
    clear tmp3
end
