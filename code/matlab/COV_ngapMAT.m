function covmat=COV_ngapMAT(smMatrix, ngap)
% 利用前面计算 得分矩阵 smMatrix，计算gapsm矩阵 (gappssm, gappcsm, gapiffr)
% 参数nagp=0时，GapIFFR 等价于 IFFR, 类似有 covpcsm, covpssm
% 输入参数：smMatrix, 跳空阶数(ngap)从0到min(len_seq)-1连续变化,
% wenjie (2022.10.04)

SM = smMatrix;
SM = 1./(1+exp(-SM));   % sigmoid归一化

% 对SM矩阵进行外滑窗(跳空n-gap)操作，如：在进化信息pssm基础上考虑部分结构信息
if ngap==0
    newSM = SM;
else
    [m,n] = size(smMatrix);
    temp = zeros(1,n);
    kk = m-ngap;        % 序列长度要大于跳空距离
    if kk<1     % 序列长度m - 跳空距离 要大于1
        error('the gap distance (ngap) is error!!')
    end
    for i=1:kk
        temp(i,:)=mean(SM([i,i+ngap],:));       
    end
    newSM=temp;
end

% 计算COV
covmat = newSM'*newSM;


