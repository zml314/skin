# skin
## background ##
lupus is a severe autoimmune disease.Clinically, LE is mainly divided into two types.One type is known as cutaneous lupus erythematosus (CLE), which mainly presents dermatological injuries and does not involve systemic damage. The other type of LE is called systemic lupus erythematosus (SLE), which involves systemic manifestations, including cutaneous, respiratory, renal, cardiovascular, and so on. Discoid lupus erythematosus (DLE) accounts for more than 80% of cases of CLE and is the most common type of CLE. On the one hand, approximately 5% of CLE patients will convert to SLE; on the other hand, approximately 20% of SLE patients show CLE lesions at the time of diagnosis or in the following years after diagnosis. Although the manifestations of cutaneous lesions are distinct between DLE and SLE, the cell composition and the underlying molecular events in cutaneous lesions of DLE and SLE remain unclear. Here, we seperated epidermis and dermis of skin lesion of DLE and SLE to find the distinction between this two types of lupus. This workflow was used for downstream analysis of single cell transcriptional.
## functions ##
Epidermis file was a workflow of single cell RNA analysis for all epidermal samples, subclustering.xml was a workflow for subclustering analysis of main cell type such as T, B, Macro/DC and NK cells, subkera.xml was used for subclustering analysis of keratinocytes.

Dermis file was a workflow of single cell RNA analysis for all epidermal samples,subfib.xml was used for subclustering analysis of fibroblasts.

## previous preparations ##
before this workflow, you may need to obtain three files after mapping rawdata to GRCh38 by cellranger, and removed doublets by Scrublet, install R and Rstudio for this analysis. Besides, several packages would be installed, such as seurat, ggplot2, matrix,dplyr, monocle, harmony, clusterprofiler, org.HS.eg.db.

## contributions ##
zhengmeiling and huzhi contributed to this workflow, if you have any questions, please contact with email in zhengmeiling314@csu.edu.cn and huzhi97@csu.edu.cn.
