% 批量生成nr下载地址，用于下载工具wget自动批量下载，
% 生成用于PSSM的库文件NR。
% wenjie (2020.04.12)

clear, clc
fidin=fopen('fuse20200412.xml','r');
nline=0;
k=0;
while ~feof(fidin)
    tline=fgetl(fidin);
    % 判断当前行是否包含(字符 'nr')，有，则该行信息抽取并写入temp.txt文档
    if  contains(tline, 'nr')
        k = k+1;
        fid1=fopen('temp.txt','a');     % 'a'--末尾追加数据方式打开文档
        
        % 用正则化，发现 URLs
        url = regexpi(tline, ...
            ['((http|https|ftp|file)://|www\.|ftp\.)',...
            '[-A-Z0-9+&@#/%=~_|$?!:,.]*[A-Z0-9+&@#/%=~_|$]'], 'match');

        fprintf(fid1, '%s\n',url{1});
        fclose(fid1);
    end
    
    nline=nline+1;
end

fclose(fidin);
