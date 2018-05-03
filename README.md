# Similarity Regression Fusion for Integrating Multi-Omics Data to Identify Cancer Subtypes

Requirements:
------------
1. python 3.6.1
   
2. numpy (1.12.1)
   
3. scikit-learn (0.18.1)

Installation:
------------
The source code can be directly called from python.


Usage:
--------------------------------
python SRF.py 

    -h|--help: get help
  
    -g|--genefile: the gene expression data file path
  
    -m|--methyfile: the methylation expression data file path
  
    -r|--mirnafile: the miRNA expression data file path
  
    -c|--clusternum: the number of clusters
  
    -x|--weight1: the weight of gene expression data
  
    -y|--weight2: the weight of methylation expression data
  
    -z|--weight3: the weight of microRNA expression data
  
       [NOTE: (x + y + z) = 1]
    
The examples of input files are available with test_net.txt, test_label.txt

Example:
--------------------------------
python SRF.py -g Gene_Expression.txt -m Methy_Expression.txt -r Mirna_Expression.txt -c 3 -x 0.5 -y 0.3 -z 0.2 

#Input:

1. Gene_Expression_Data.txt

2. Methy_Expression_Data.txt

3. Mirna_Expression_Data.txt

#Output:
Patient Subtypes Labels. The first column represents patients' ID, and the second column represents patients' label.

      
   
----------------------------------
Copyright and License Information
----------------------------------
Copyright (C) 2018 Northwestern Polytechnical University, Xi’an, China. Yang Guo(gyang@mail.nwpu.edu.cn), Jianning Zheng(jennings@mail.nwpu.edu.cn)
