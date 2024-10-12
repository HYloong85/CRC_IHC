#### CRC_IHC_staining_identification.m 

It is a demo of staining identification visualization and staining ratio calculation. Running this script can obtain visualized results of staining identification (**TMAdemo_stain.tiff** ) and the proportion of staining in various tissues (**tissuePercent_TMAdemo.mat**).

#### nor_and_tum.m

It is a script that distinguishes between normal stroma and tumor stroma.

#### net.mat

It is the softmax model weight for staining identification.

**outout_masks** and **outout_margin** are the visualization results and edge information saved by tissue classification.

**mask_new**  save the mask of five types of tissues classification results (glands, normal stroma, tumor stroma, tumors, and others) with added normal stroma and tumor stroma using a mask

**img_new_mask_out_TMAdemo.tiff**  Visualization of classification results for five types of organizations