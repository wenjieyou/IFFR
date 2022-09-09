% 把seq_FASTA目录下单条序列，合并成多序列的FASTA格式文件，
% 用于website中计算profeat或者PseAAC等操作
% wenjie (2016.06.20)

clear, clc

fid1=fopen('temp222.txt', 'a');     % 'a'--末尾追加数据方式打开文档

for k=1:2306
    fn=['alter_seq_FASTA\seq_FASTA', num2str(k), '.fasta'];
    fidin=fopen(fn, 'r');
    nline=0;    
   
    while ~feof(fidin)
        tline=fgetl(fidin);       
        fprintf(fid1, '%s\n',tline);        
        nline=nline+1;
    end
    
    fclose(fidin);
end
fclose(fid1);

