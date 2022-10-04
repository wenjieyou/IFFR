% 特征表示：GapIFFR
% wenjie (2022.10.04)

tic,
clear, clc
load PDB1075_p2sm26Matrix
numseq = size(p2sm,1);
res=zeros(numseq,1);

for k=1:numseq
    lenseq(k,:) = size(p2sm{k},1);
end
minlen = min(lenseq);

for i=0:minlen-1
    parfor k=1:numseq
        Matrix=p2sm{k};
        temp=COV_ngapMAT(Matrix,i);
        cov_vec(k,:)=temp(:)';
    end
    res = [res, cov_vec];
end
res(:,1) = [];
cov_ngapp2sm=res;

numgap = minlen;    % gap参数：gap=0无跳空操作；gap=1相邻残基；gap>1不相邻残基；
%% 以下对角化操作
aa=ones(size(p2sm{1},2), size(p2sm{1},2));
bb=triu(aa);
cc=bb(:)';
cc=logical(cc);
cov_ngapp2sm=cov_ngapp2sm(:,repmat(cc,1,numgap));

toc
filename = ['PDB1075_vec351_gapp2sm26_gap0_', num2str(minlen-1), '.mat'];
save(filename, 'cov_ngapp2sm');
