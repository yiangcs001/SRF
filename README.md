# Similarity Regression Fusion for Integrating Multi-Omics Data to Identify Cancer Subtypes

Requirements:
------------
   1.python 3.6.1
   
   2.numpy (1.12.1)
   
   3.scikit-learn (0.18.1)

Installation:
------------
  The source code can be directly called from python.


Usage:
--------------------------------
  -h              : Get Help.
  
  -g/--genefile   : The Gene Expression Data File Path.
  
  -m/--methyfile  : The Methylation Expression Data File Path.
  
  -r/--mirnafile  : The microRNA Expression Data File Path.
  
  -c/--clusternum : The Number of Clusters.
  
  -x/--weight1    : The Weight of Gene Expression Data.
  
  -y/--weight2    : The Weight of Methylation Expression Data.
  
  -z/--weight3    : The Weight of microRNA Expression Data.
  
  [NOTE: (x + y + z) = 1]
    
The examples of input files are available with test_net.txt, test_label.txt. 

Example:
--------------------------------
python SRF.py -g GBM_Gene_Expression.txt -m GBM_Methy_Expression.txt -r GBM_Mirna_Expression.txt -c 3 -x 0.5 -y 0.3 -z 0.2

#Input:

1.Gene_Expression_Data.txt

2.Methy_Expression_Data.txt

3.Mirna_Expression_Data.txt


#Output:
1.Patient Subtypes Labels.The first column represents patients' ID, and the second column represents patients' label.

      
   
----------------------------------
Copyright and License Information
----------------------------------
Copyright (C) 2017 Northwestern Polytechnical University, Xi’an, China. Yang Guo(gyang@mail.nwpu.edu.cn), Jianning Zheng(jennings@mail.nwpu.edu.cn)
