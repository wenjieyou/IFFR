% ��������nr���ص�ַ���������ع���wget�Զ��������أ�
% ��������PSSM�Ŀ��ļ�NR��
% wenjie (2020.04.12)

clear, clc
fidin=fopen('fuse20200412.xml','r');
nline=0;
k=0;
while ~feof(fidin)
    tline=fgetl(fidin);
    % �жϵ�ǰ���Ƿ����(�ַ� 'nr')���У��������Ϣ��ȡ��д��temp.txt�ĵ�
    if  contains(tline, 'nr')
        k = k+1;
        fid1=fopen('temp.txt','a');     % 'a'--ĩβ׷�����ݷ�ʽ���ĵ�
        
        % �����򻯣����� URLs
        url = regexpi(tline, ...
            ['((http|https|ftp|file)://|www\.|ftp\.)',...
            '[-A-Z0-9+&@#/%=~_|$?!:,.]*[A-Z0-9+&@#/%=~_|$]'], 'match');

        fprintf(fid1, '%s\n',url{1});
        fclose(fid1);
    end
    
    nline=nline+1;
end

fclose(fidin);
