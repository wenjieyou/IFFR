function covmat=COV_ngapMAT(smMatrix, ngap)
% ����ǰ����� �÷־��� smMatrix������gapsm���� (gappssm, gappcsm, gapiffr)
% ����nagp=0ʱ��GapIFFR �ȼ��� IFFR, ������ covpcsm, covpssm
% ���������smMatrix, ���ս���(ngap)��0��min(len_seq)-1�����仯,
% wenjie (2022.10.04)

SM = smMatrix;
SM = 1./(1+exp(-SM));   % sigmoid��һ��

% ��SM��������⻬��(����n-gap)�������磺�ڽ�����Ϣpssm�����Ͽ��ǲ��ֽṹ��Ϣ
if ngap==0
    newSM = SM;
else
    [m,n] = size(smMatrix);
    temp = zeros(1,n);
    kk = m-ngap;        % ���г���Ҫ�������վ���
    if kk<1     % ���г���m - ���վ��� Ҫ����1
        error('the gap distance (ngap) is error!!')
    end
    for i=1:kk
        temp(i,:)=mean(SM([i,i+ngap],:));       
    end
    newSM=temp;
end

% ����COV
covmat = newSM'*newSM;


