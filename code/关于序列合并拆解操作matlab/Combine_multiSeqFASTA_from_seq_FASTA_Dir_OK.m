% ��seq_FASTAĿ¼�µ������У��ϲ��ɶ����е�FASTA��ʽ�ļ���
% ����website�м���profeat����PseAAC�Ȳ���
% wenjie (2016.06.20)

clear, clc

fid1=fopen('temp222.txt', 'a');     % 'a'--ĩβ׷�����ݷ�ʽ���ĵ�

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

