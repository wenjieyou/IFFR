% �Ѷ�����fasta�ļ��ָ�ɵ������ļ���ʽ��
% ��ԭʼD8244.xls�е�PDB_code����д��fasta�����ļ���
% wenjie (2019.07.30)

clear, clc
fidin1=fopen('temp.txt','r');
fidin2=fopen('D8244.csv','r');
nline=0;
k=0;
while ~feof(fidin1)
    tline=fgetl(fidin1);
    
    % �жϵ�ǰ���Ƿ�Ϊ������غ�(���ַ� '>' Ϊ���)
    if tline(1)=='>'
        k = k+1;
        fn=['temp\seq_FASTA', num2str(k), '.fasta'];		% �ڵ�ǰĿ¼���½���Ŀ¼temp���ָ��ĵ��������д���ڴˡ�
        fid1=fopen(fn,'a');     % 'a'--ĩβ׷�����ݷ�ʽ���ĵ�
        tline_title=fgetl(fidin2);
        % ��غ�( '>' )��д��PDB�ͽṹ����Ϣ
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
