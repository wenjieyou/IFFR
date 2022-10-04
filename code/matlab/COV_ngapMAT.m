function covmat=COV_ngapMAT(pssmMatrix, ngap)
% 利用前面计算psiblast得到的pssmMatrix矩阵，计算gappssm矩阵
% 输入参数：pssmMatrix, 跳空阶数(ngap)从0到min(len_seq)-1连续变化,
% wenjie (2022.10.04)

blastPSSM = pssmMatrix;
blastPSSM = 1./(1+exp(-blastPSSM));   % sigmoid归一化

% 对PSSM矩阵进行外滑窗(跳空n-gap)操作，在进化信息pssm基础上考虑部分结构信息
if ngap==0
    newPSSM=blastPSSM;
else
    [m,n] = size(pssmMatrix);
    temp = zeros(1,n);
    kk = m-ngap;        % 序列长度要大于跳空距离
    if kk<1     % 序列长度m - 跳空距离 要大于1
        error('the gap distance (ngap) is error!!')
    end
    for i=1:kk
        temp(i,:)=mean(blastPSSM([i,i+ngap],:));       
    end
    newPSSM=temp;
end

% 计算COV
covmat = newPSSM'*newPSSM;


