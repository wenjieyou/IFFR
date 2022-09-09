% Ԭ���ṩ�ĵ����ʶ�����ʽ���ݣ�����PSSM
% ��python���ݸ�ʽ.npy�������ļ���ת��Ϊmatlab���ݸ�ʽ.mat
%       from scipy import io
%       mat=np.load('CullPDB13114.npy')
%       io.savemat('seq.mat',{'seq':mat})
% Ȼ��matlab����Щ����д��fasta��ʽ
% wenjie (2081.11.24)

clear, clc
load seq.mat

for i=1:size(seq,1)
    tmp1 = seq(i,:,:);          % ��ά��������ʽ
    tmp2 = shiftdim(tmp1);      % ��ά�г�ȡһ����ά��
    tmp3 = tmp2(1,:);               % ��ȡ����һά����ԭʼ����
    tmp3(find(isspace(tmp3))) = [];     % ɾ���м���ܵĿո�
    len = length(tmp3);
    
    filename = ['seq_FASTA\seq_FASTA_', num2str(i)];
    fid = fopen(filename,'a');
    firstline = ['>', num2str(i)];          % ����д��fasta��ʽ��'>'
    fprintf(fid, '%s\n', firstline);
    
    int = fix(len/60);          %�����г���ȡ����
    rem = mod(len, 60);     %�����г�������
    if int >1
        for j = 1:fix(len/60)
            tline = tmp3((j-1)*60+1:j*60);
            fprintf(fid, '%s\n', tline);
        end
    end
    if rem > 0
        fprintf(fid, '%s', tmp3(end+1-rem:end));
    end
    fclose(fid);
    clear tmp3
end
