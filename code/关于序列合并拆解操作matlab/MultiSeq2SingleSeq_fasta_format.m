% �Ѷ�����fasta�ļ��ָ�ɵ������ļ���ʽ�����ڼ���pssm����
% wenjie (2017.07.17)

clear, clc
fidin=fopen('ZW225_seq_FASTA.seq','r');
nline=0;
k=0;
while ~feof(fidin)
    tline=fgetl(fidin);
    
    % �жϵ�ǰ���Ƿ�Ϊ������غ�(���ַ� '>' Ϊ���)
    if tline(1)=='>'
        k = k+1;
        fn=['temp\seq_FASTA', num2str(k), '.fasta'];		% �ڵ�ǰĿ¼���½���Ŀ¼temp���ָ��ĵ��������д���ڴˡ�
        fid1=fopen(fn,'a');     % 'a'--ĩβ׷�����ݷ�ʽ���ĵ�
        fprintf(fid1, '%s\n',tline);
        fclose(fid1);
    else
        if length(tline)>60
            nn=floor(length(tline)/60);
            fid2=fopen(fn,'a');     % 'a'--ĩβ׷�����ݷ�ʽ���ĵ�
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
