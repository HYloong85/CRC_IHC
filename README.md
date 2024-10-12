# CRC-IHC

Colorectal cancer (CRC) is one of the most common gastrointestinal malignancies, and the prognosis of advanced patients is usually poor. This project collected tissue microarray pathological images of 154 postoperative patients with colorectal cancer and studied the impact of infiltration expression of 15 immune markers on prognosis.
This project establishes an immune automatic scoring system applicable to different tissue types, classifies pathological tissues through deep learning networks, and calculates the percentage of stained pixels as the immune score. In survival analysis, overall survival (OS) and disease-free survival (RFS) were used as endpoints for this study, and X-tile was used to select cutoff values for dividing patient cohorts.

## Experiment 1: Tissue classification

For tissue classification, we trained a patch-based CNN model. This model categorized tissues into four types: gland, tumor, stroma, and all other types (including muscle layer, debris, and background).  

The folder "Tissue Classification" contains a demo of tissue classification, which can be used to obtain the proportion of four types of tissues in the image and visualize the classification results using code (CRC_IHC_tissue_classification.m) and model weights(VGG19_model.mat).



## Experiment 2: Staining identification

For staining identification, we trained pixel-based machine learning classification models to categorize pixels into three classes: stained tissue, unstained tissue, and background. 

The folder "Staining Identification" contains a demo of staining identification. You can use code (CRC_IHC_staining_identification.m) and model weights (net.mat) to obtain the proportion of staining in each tissue in the image and visualize the results of staining classification.



## Experiment 3: Prognostic analysis of individual immune markers

For each of the 15 markers (CD3, CD4, CD8, CD20, CD45RO, Granzyme B, S100, CD68, Tryptase, CD57, FOXP3, HLA-DR, Fas, FasL, and IL-17), we calculated the proportion of stained pixels in eight tissue types as the IHC scores. These tissue types include glands and stroma from both healthy and paracancerous tissues, as well as tumor and stroma from both invasive margin and tumor center. The combination of 8 tissue types and 15 immune markers yielded a total of 120 IHC scores for survival analysis. We used the log-rank test to determine the significance of survival differences between patients in the high and low expression groups.

The folder "Prognostic analysis of individual immune markers" contains IHC scores for 15 immune markers from 8 different types of tissues (completed_maker *.csv), the cutoff value obtained from xtile (OS/RFS xtile cutoff.csv) and prognostic analysis code(univariate logranktest xtile.R). 



## Experiment 4:  Prognostic analysis of combined marker score

We calculated a combined marker score based on their expression patterns.  For each survival-associated marker, one point was added to the combined score for patients exhibiting favorable expression patterns: high expression of a positive marker or low expression of a negative marker. All other patterns were assigned zero points. The total points were summed, and patients were stratified into high and low groups based on the median of the combined score.

The folder "Prognostic analysis of combined marker score" contains IHC scores for 15 immune markers from 8 different types of tissues(completed_maker *.csv), the cutoff value obtained from xtile(completed_maker *.csv) and code for calculating the combined marker score(combined score.R). KM survival curve can be obtained through prognostic analysis code.



## Experiment 5:  Lasso Cox of THIR score

We developed a Tumor-to-Healthy Immune Ratio (THIR) to compare immune marker expression in the stroma of tumor versus healthy tissues. The THIR score was calculated to quantify the differential immune infiltration between the tumor and healthy regions: 
$$
THIR=(IHC score in the stroma of invasive margin+IHC score in the stroma of tumor center)/(IHC score in the stroma of normal tissue+IHC score in the stroma of paracancerous tissue)
$$


The folder "Lasso Cox of THIR score" contains data from 8 types of tissues combined (all fea data.csv) and lasso cox analysis code (lasso cox.R). KM survival curve can be obtained through lasso cox analysis code. The timeROC curve can be plotted through code (timeROC.R).



## How to run the code?

The code is written in MATLAB and R. For the five abovementioned experiments, run the scripts in the corresponding folders sequentially according to the file names. The folder "tools"  contains the tool packages used.



## Contact information

If you have any questions, feel free to contact me.

Yulong Han

Shenzhen University, Shenzhen, China
E-mail: hanyl0805@163.com