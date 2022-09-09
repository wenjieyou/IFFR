function replaceX(FileName)
% 数据预处理：复制序列，把原序列中的异常碱基用随机生成的碱基替换
% NN 替换掉的数目
% wenjie modified (2015.03.21)

prt='NDRIWFCAPSHQEKLYGTVM';  % pre
fpr=fopen(FileName,'rt'); % 使用 CD-HIT 处理过的序列文件
fpw=fopen(strcat('tmp_',FileName),'wt'); % OutPut File
%%
NN=0;
while ~feof(fpr)
    tline=fgetl(fpr);
    if isempty(tline) || tline(1)=='>'      % 空行和名称行直接复制
        fprintf(fpw,'%s\n',tline);
        continue;
    end
    % 非空行和名称行
    n=length(tline);
    for i=1:n
        if isempty(find(tline(i)==prt, 1)) % 若此字符是异常字符
            nnx=ceil( 20*rand(1) );
            tline(i)=prt(nnx);  % 以等概率选取一个氨基酸代替
            NN=NN+1;
        end
    end
    fprintf(fpw,'%s\n',tline);
end
fclose(fpr);  
fclose(fpw);
fprintf('Successfully! %d chars to be replaced!\n', NN);

