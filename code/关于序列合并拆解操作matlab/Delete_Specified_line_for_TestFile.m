% ��website���ɵ�PseAAC�ļ���ɾ�����з���ֵ���ַ�������PseAAC420��ֵ��
% ÿ21��Ϊһ�����е���Ϣ�����ҵ�һ��ΪAAC20
% wenjie (2016.06.20)

clear, clc
fidin=fopen('Indset2_binder_PseAAC420_from_website.txt','r');
nline=0;
k=0;
while ~feof(fidin)
    tline=fgetl(fidin);
    
    % �жϵ�ǰ��Ϊ��������غ�(���ַ� '>' Ϊ���)�Ϳ��У���������Ϣд��temp.txt�ĵ�
    if  isempty(strfind(tline, '>'))&&length(tline)>1
        k = k+1;
        fid1=fopen('temp.txt','a');     % 'a'--ĩβ׷�����ݷ�ʽ���ĵ�
        fprintf(fid1, '%s\n',tline);
        fclose(fid1);
    end
    
    nline=nline+1;
end

fclose(fidin);
