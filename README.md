# shinyDeepGxP
Overview
This repository contains the complete shinyDeepGxP application, including source code, pretrained DeepGxP model files, reference datasets, annotation resources, interpretation results, and an example dataset required to reproduce the functionality of the public web server.
shinyDeepGxP is a user-friendly R Shiny web application for predicting surface protein abundance from single-cell RNA sequencing (scRNA-seq) data using the DeepGxP deep learning framework. The application enables researchers to infer the abundance of 224 surface proteins from transcriptomic profiles and to investigate the biological mechanisms underlying model predictions.
The application contains two complementary modules:
Explore Model
    Investigate genes contributing to the prediction of a selected surface protein and explore their associated biological pathways and functional enrichments.
Predict Proteins
    Predict the abundance of 224 surface proteins from user-uploaded scRNA-seq datasets and perform downstream analyses including clustering, cell type annotation, and UMAP visualization.


Web Server:
https://shiny.crc.pitt.edu/deepgxp/


Recommended Use
The DeepGxP model was trained using human PBMC CITE-seq data comprising 53,364 cells, 19,738 genes, and 224 surface proteins. Accordingly, predictions are expected to perform best on human blood and immune-cell datasets. Results from unrelated tissues or cell types should be interpreted with caution.
The repository contains the same pretrained models, reference datasets, and interpretation resources used by the public web server.


Running Locally
Option 1: RStudio
1. Download and extract the source code.
2. Open the project in RStudio.
3. Open either ui.R or server.R.
4. Click 'Run App' button.

Option 2: Docker
1. Download and extract the source code.
2. Install Docker.
3. Navigate to the project directory using cmd or linux backend.
4. Build the Docker image: 
		docker build --platform=linux/amd64 -t deepgxp .
5. Run the container: 
		docker run --rm -it --platform=linux/amd64 -p 7891:3838 deepgxp
6. Open the application in a web browser:
		http://localhost:7891.


Repository Contents
server.R
    Shiny server logic
ui.R
    Shiny user interface
index.html
    Web page layout
function/
    Supporting R functions
www/
    HTML help pages and web assets
data/
    Model metadata, reference datasets, annotation resources, and interpretation results
Predictions/
    Pretrained DeepGxP model files
Predictions/samples/
    Example datasets for demonstration and testing
Dockerfile
    Docker deployment configuration
README.md
    Documentation


Software Requirements
The application was tested using:
R 4.4.2

Major dependencies include:
Seurat
SeuratObject
SeuratDisk
reticulate
Rmagic
keras
tensorflow
shiny


Data Availability
All files required to run shinyDeepGxP locally, including pretrained models, reference datasets, annotation resources, interpretation results, and example datasets, are provided within this repository.


Citation:
shinyDeepGxP web server
Tsai HM, Hsiao TH, Hsu YC, Wang LJ, Chiu YC, Chuang EY, Chen Y.
shinyDeepGxP: A user-friendly R Shiny app for predicting surface protein abundance from scRNA-seq expression using deep learning in blood cells.
Manuscript under review.

DeepGxP model
Tsai HM, Hsiao TH, Chiu YC, Huang Y, Chuang EY, Chen Y.
Predicting and interpreting protein and phosphoprotein abundance from pan-cancer and single-cell transcriptomes.
iScience, 2026.


Contact information:
Yidong Chen — Cheny8@uthscsa.edu
Eric Y. Chuang — chuangey@ntu.edu.tw
Yu-Chiao Chiu — YUC250@pitt.edu


License
This software is distributed for academic and research use.