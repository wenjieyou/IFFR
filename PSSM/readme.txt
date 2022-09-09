PSI-BLAST的特色是每次用profile搜索数据库后再利用搜索的结果重新构建profile，然后用新的profile再次搜索数据库，如此反复直至没有新的结果产生为止。PSI-BLAST先用带空位的BLAST搜索数据库，将获得的序列通过多序列比对来构建第一个profile。PSI-BLAST自然地拓展了BLAST方法，能寻找蛋白质序列中的隐含模式，有研究表明这种方法可以有效的找到很多序列差异较大而结构功能相似的相关蛋白，甚至可以与一些结构比对方法，如threading相媲美。PSI-BLAST服务可以在NCBI的BLAST主页上找到，还可以从NCBI的FTP服务器上下载PSI-BLAST的独立程序。

===============================

如何用psi-blast生成pssm矩阵?
--现在想先用f[x]=1/1+exp(-x)做正规化，
再用滑动窗口法为每个残基提取w*20维特征向量，
当残基在序列两端时用0补齐。
================================

blast怎么选择其他替换矩阵
在替换矩阵中,每个位置的打分是在相关蛋白局部比对模块中观察到的替换的频率
===============================
一个完整的psi-blast命令如下

psiblast -query 1_dp.fasta -db nr -out 1.txt -num_iterations 3 -out_ascii_pssm 1.pssm

C:\Program Files\NCBI\blast-2.2.30+\bin\psiblast -db C:\\Blast\\db\\nr.00 -in_msa seq_FASTA1.fasta -evalue 0.001 -num_iterations 3 -out PSSM1.txt -out_ascii_pssm PSSMmatrix1.txt

psiblast -in_msa seq_FASTA2.fasta -out_ascii_pssm pssm2.txt -subject 1

psiblast -db nr.00 -in_msa seq_FASTA2.fasta -evalue 0.001 -num_iterations 3 -out PSSM2.txt -out_ascii_pssm PSSMmatrix2.txt

 *** PSSM engine options
 -in_msa <File_In>
   File name of multiple sequence alignment to restart PSI-BLAST
    * Incompatible with:  in_pssm, query, query_loc, phi_pattern

-query <File_In>
   Input file name
   Default = `-'
    * Incompatible with:  in_msa, msa_master_idx, ignore_msa_master, in_pssm

-out_pssm <File_Out>
   File name to store checkpoint file
 -out_ascii_pssm <File_Out>
   File name to store ASCII version of PSSM

===========20200418=================================

PSI-BLAST(Position-Specific Iterative Basic Local Alignment Search Tool)是blastp的一种，用于发现与查询序列更远距离相关的序列。
PSSM（position-specific scoring matrix），PSI-BLAST中要用到的中间结果，也被称为profile（配置文件）。PSSM捕获保持模式的对齐并将其存储为对齐中每个位置的分数矩阵 - 高度保守的位置接收高分数，弱保守位置接收分数接近零。该配置文件用于代替原始替换矩阵，以进一步搜索数据库以检测与PSSM指定的保护模式匹配的序列。

步骤
PSI-BLAST一般分为以下几个步骤：

1.用提交的序列在数据库中进行搜索，（与blastp相同）；
2.对初始的结果进行多序列比对，构建初始profile，即PSSM；
3.用PSSM作为query，搜索数据库，发现较远的相似性序列，并给出E-value；
4.构建新的PSSM，重复步骤3。这个过程不断迭代，直到没有符合要求(E-value threshold,比如0.005,有时可放宽至0.01)的新序列或者到达迭代次数要求。

运行
../ncbi-blast-2.2.25+/bin/psiblast  -query ../protein/1A01_HUMAN.fasta -evalue .001 -inclusion_ethresh .002 -db ../databases/nr90-2012/nr90 -num_iterations 3  -seg yes -outfmt '7 std qseq sseq stitle' -out 1A01_HUMAN_100.output -max_target_seqs 100
../ncbi-blast-2.2.25+/bin/psiblast  -query ../protein/1A01_HUMAN.fasta -evalue .001 -inclusion_ethresh .002 -db ../databases/nr90-2012/nr90 -num_iterations 3  -seg yes -outfmt '7 std qseq sseq stitle' -out 1A01_HUMAN_100.output -max_target_seqs 500
-query:搜索的序列；
-num_iterations：迭代次数；
-seg:是否SEG过滤；
-max_target_seqs：最大匹配数量，用100、500做了两次
-outfmt:输出格式设置，改了一下方便后续处理。
---------------------------------------------------------

PSSM中同一个氨基酸在不同位点的得分不同，因此比PAM和BLOSUM更强调保守性和结构性，而不是序列本身.


