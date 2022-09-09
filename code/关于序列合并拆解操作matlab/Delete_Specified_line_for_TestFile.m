% 把website生成的PseAAC文件，删除所有非数值型字符，保留PseAAC420数值，
% 每21行为一个序列的信息，并且第一行为AAC20
% wenjie (2016.06.20)

clear, clc
fidin=fopen('Indset2_binder_PseAAC420_from_website.txt','r');
nline=0;
k=0;
while ~feof(fidin)
    tline=fgetl(fidin);
    
    % 判断当前行为非序列入藏号(以字符 '>' 为标记)和空行，将该行信息写入temp.txt文档
    if  isempty(strfind(tline, '>'))&&length(tline)>1
        k = k+1;
        fid1=fopen('temp.txt','a');     % 'a'--末尾追加数据方式打开文档
        fprintf(fid1, '%s\n',tline);
        fclose(fid1);
    end
    
    nline=nline+1;
end

fclose(fidin);
