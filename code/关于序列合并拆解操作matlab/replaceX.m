function replaceX(FileName)
% ����Ԥ�����������У���ԭ�����е��쳣�����������ɵļ���滻
% NN �滻������Ŀ
% wenjie modified (2015.03.21)

prt='NDRIWFCAPSHQEKLYGTVM';  % pre
fpr=fopen(FileName,'rt'); % ʹ�� CD-HIT ������������ļ�
fpw=fopen(strcat('tmp_',FileName),'wt'); % OutPut File
%%
NN=0;
while ~feof(fpr)
    tline=fgetl(fpr);
    if isempty(tline) || tline(1)=='>'      % ���к�������ֱ�Ӹ���
        fprintf(fpw,'%s\n',tline);
        continue;
    end
    % �ǿ��к�������
    n=length(tline);
    for i=1:n
        if isempty(find(tline(i)==prt, 1)) % �����ַ����쳣�ַ�
            nnx=ceil( 20*rand(1) );
            tline(i)=prt(nnx);  % �Եȸ���ѡȡһ�����������
            NN=NN+1;
        end
    end
    fprintf(fpw,'%s\n',tline);
end
fclose(fpr);  
fclose(fpw);
fprintf('Successfully! %d chars to be replaced!\n', NN);

