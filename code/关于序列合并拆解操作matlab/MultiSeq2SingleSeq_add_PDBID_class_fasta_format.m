% 把多序列fasta文件分割成单序列文件形式，
% 对原始D8244.xls中的PDB_code数据写入fasta数据文件中
% wenjie (2019.07.30)

clear, clc
fidin1=fopen('temp.txt','r');
fidin2=fopen('D8244.csv','r');
nline=0;
k=0;
while ~feof(fidin1)
    tline=fgetl(fidin1);
    
    % 判断当前行是否为序列入藏号(以字符 '>' 为标记)
    if tline(1)=='>'
        k = k+1;
        fn=['temp\seq_FASTA', num2str(k), '.fasta'];		% 在当前目录下新建子目录temp，分割后的单条子序列存放于此。
        fid1=fopen(fn,'a');     % 'a'--末尾追加数据方式打开文档
        tline_title=fgetl(fidin2);
        % 入藏号( '>' )后写入PDB和结构类信息
        fprintf(fid1, '%s\n',[tline, ': [PDB_Code,Class,Region]=', tline_title]);   
        fclose(fid1);
    else
        fid3=fopen(fn,'a');
        fprintf(fid3, '%s\n',tline);
        fclose(fid3);
    end
    
    nline=nline+1;
end

fclose(fidin1);
fclose(fidin2);
