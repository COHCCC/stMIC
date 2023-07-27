# stMIC: Spatial Transcriptomics and Morphological Integrated Clustering

Spatial transcriptomics (ST) represents a pivotal advancement in biomedical research, enabling the transcriptional profiling of cells within their morphological context and providing a pivotal tool for understanding spatial heterogeneity in cancer tissues. However, current analytical approaches, akin to single-cell analysis, largely depend on gene expression, underutilizing the rich morphological information inherent in the tissue. 

We present a novel method integrating spatial transcriptomics and histopathological image data to better capture biologically meaningful patterns in patient data, focusing on aggressive cancer types such as glioblastoma and triple-negative breast cancer. We used a ResNet-based deep learning model to extract key morphological features from high-resolution whole-slide histology images. We merged these with the gene expression matrix to create a comprehensive dataset, followed by unsupervised clustering. Our integrative approach successfully pinpointed key biological features identified by manual histopathologies, such as regions of fibrosis and necrosis regions, as well as improved edge definition in EGFR-rich areas. Importantly, our combinatorial approach unveiled new features evident in histopathology, which were missed by gene-expression-only analysis.

<img src="/Users/ninasong/Desktop/spatialProject/PBS2024/stMIC/plot/workflow.png">



For the step-by-step tutorial with explanation, please refer to: [*run example*](https://github.com/USCDTG/stMIC/blob/main/Resnet_gene_exp_integration.ipynb)

For result visulization, please refer to 



