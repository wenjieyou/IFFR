function covmat=COV_ngapMAT(pssmMatrix, ngap)
% ����ǰ�����psiblast�õ���pssmMatrix���󣬼���gappssm����
% ���������pssmMatrix, ���ս���(ngap)��0��min(len_seq)-1�����仯,
% wenjie (2022.10.04)

blastPSSM = pssmMatrix;
blastPSSM = 1./(1+exp(-blastPSSM));   % sigmoid��һ��

% ��PSSM��������⻬��(����n-gap)�������ڽ�����Ϣpssm�����Ͽ��ǲ��ֽṹ��Ϣ
if ngap==0
    newPSSM=blastPSSM;
else
    [m,n] = size(pssmMatrix);
    temp = zeros(1,n);
    kk = m-ngap;        % ���г���Ҫ�������վ���
    if kk<1     % ���г���m - ���վ��� Ҫ����1
        error('the gap distance (ngap) is error!!')
    end
    for i=1:kk
        temp(i,:)=mean(blastPSSM([i,i+ngap],:));       
    end
    newPSSM=temp;
end

% ����COV
covmat = newPSSM'*newPSSM;


