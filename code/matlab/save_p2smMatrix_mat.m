
clear, clc

load PDB1075_uniref50_pssmMatrix
load PDB1075_LiuBin3_pcsm6Matrix
seqnum=size(pssm,1);

for k=1:seqnum
    tmp1=pssm{k};
    tmp2=pcsm{k};
    res{k,:}=[tmp1,tmp2];
end

p2sm=res;
save PDB1075_p2sm26Matrix p2sm label