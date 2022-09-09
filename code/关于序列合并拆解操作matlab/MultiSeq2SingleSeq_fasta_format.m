% 把多序列fasta文件分割成单序列文件形式，用于计算pssm矩阵
% wenjie (2017.07.17)

clear, clc
fidin=fopen('ZW225_seq_FASTA.seq','r');
nline=0;
k=0;
while ~feof(fidin)
    tline=fgetl(fidin);
    
    % 判断当前行是否为序列入藏号(以字符 '>' 为标记)
    if tline(1)=='>'
        k = k+1;
        fn=['temp\seq_FASTA', num2str(k), '.fasta'];		% 在当前目录下新建子目录temp，分割后的单条子序列存放于此。
        fid1=fopen(fn,'a');     % 'a'--末尾追加数据方式打开文档
        fprintf(fid1, '%s\n',tline);
        fclose(fid1);
    else
        if length(tline)>60
            nn=floor(length(tline)/60);
            fid2=fopen(fn,'a');     % 'a'--末尾追加数据方式打开文档
            for i=1:nn
                fprintf(fid2, '%s\n',tline((i-1)*60+1:(i-1)*60+60));
            end
            fprintf(fid2, '%s\n',tline(nn*60+1:end));
            fclose(fid2);
        else
            fid3=fopen(fn,'a');
            fprintf(fid3, '%s\n',tline);
            fclose(fid3);
        end
    end
    
    nline=nline+1;
end

fclose(fidin);
