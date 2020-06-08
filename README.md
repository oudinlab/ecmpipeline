# ecmpipeline
This repository includes the Matlab code and pipeline data used in the Oudin Lab. Please see Baskaran et al., APL Bioengineering, 2020 for more information.

Files Included:
- Matlab Code Notes - word document with annotations on basic code elements and how to use the scripts (also links to resources to learn more about the analyses)

- PLSData.xlsx: data for PLSR (averages for each variable)

Matlab Scripts and Functions:
- PLSRCode.m: performs PLSR analysis (outputs R2, loadings and scores plots, VIP Scores, and predictions)
- PLS_CV.m: performs a cross validation on PLS model and calculates Q2 metric
- CrossValFunc.m: performs the leave-one-out cross validation (called by PLS_CV.m and PermutationTest.m)
- q2calc.m: calculates Q2 value (called by PLS_CV.m and PermutationTest.m)
