# RMA tracing
This repository contains Python code to process the output data from Compound Discoverer to indentify compounds that has incorporated an isotopically labelled substrate.

For example, mammalian cell lines grown in media with 50 % unlabelled Cysteine and 50 % 13C3, 15N1 Cysteine (m+4) will incorporate the same labelling fraction into downstream metabolites like glutathione and lactoylglutathione.


## Examples
Three examples are provided:
1. Cysteine and glutamine tracing in three cancer cell lines. Input and output in the folder `projects/three-cell-lines` and notebook to process the data: `peak-pair-analysis_three-cell-lines.ipynb`
2. Cysteine tracing in two cancer cell lines treated with six alkylating agents. Input and output in the folder `projects/alkylation-agents` and notebook to process the data: `peak-pair-analysis_alkylation-agents.ipynb`
3. Cysteine tracing in a set of bile duct cancer cell lines. Input and output in the folder `projects/bile-duct-cells_cys-tracing` and notebook to process the data: `peak-pair-analysis_bile-duct-cells_cys-tracing.ipynb`


## Input requirements
* Excel file with the output of Compound Discoverer. This file should contain one row per compound and three columns for: "Name", "Molecular Weight", and "RT \[min\]". And then followed by the columns for the peak area per sample e.g.: "Area: KD031720_031720_BD_C01.raw (F1)", "Area: KD031720_031720_BD_C02.raw (F2)", etc. The peak area columns must contain the tag "Area", which it will automatically if the Compound Discoverer output is not altered.
* A JSON file describing each input sample, its type, label etc. See the examples provided in the `projects` folder.






