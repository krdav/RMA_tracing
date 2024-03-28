# RMA tracing
This repository contains Python code to process the output data from Compound Discoverer to indentify compounds that has incorporated an isotopically labelled substrate.

For example, mammalian cell lines grown in media with 50 % unlabelled Cysteine and 50 % 13C3, 15N1 Cysteine (m+4) will incorporate the same labelling fraction into downstream metabolites like glutathione and lactoylglutathione.


## Examples
Three examples are provided:
1. Cysteine and glutamine tracing in three cancer cell lines. Input and output in the folder `projects/three-cell-lines` and notebook to process the data: `peak-pair-analysis_three-cell-lines.ipynb`
2. Cysteine tracing in two cancer cell lines treated with six alkylating agents. Input and output in the folder `projects/alkylation-agents` and notebook to process the data: `peak-pair-analysis_alkylation-agents.ipynb`
3. Cysteine tracing in a set of bile duct cancer cell lines. Input and output in the folder `projects/bile-duct-cells_cys-tracing` and notebook to process the data: `peak-pair-analysis_bile-duct-cells_cys-tracing.ipynb`


## How to run
1. Clone this repository.
2. Install Python packages that are dependencies (see below).
3. Run one or several of the example notebooks.
4. If you want to use this code with your own data start by copying one of the example notebooks and use this as a boilerplate.
5. Prepare your input data according to the requirements described below. We recommend copying the folder structure of one of the examples e.g. `projects/bile-duct-cells_cys-tracing` and use these input files as templates.
6. Adjust the parameters for peak filtering and peak pair finding. These are input dependent e.g. the mass shift and retention time tolerances must reflect the mass accuracy and chromatographic setup.
7. Adjust the use of tracer (see below) if different from those used in the examples.
8. Adjust the peak pair area ratio to reflect the mixing ratio of labelled and unlabelled tracer. This can also be performed empirically by inputting expected peak pairs to the `pick_ratio` method shown in the example notebooks.
9. Generate peak pair using the `find_pairs` method as shown in the example notebooks.
10. Inspect peak pairs and adjust filtering parameters. Repeat filter adjustments until satisfied.


### Choice of tracer
We provide examples using two tracers: m+4 cysteine (13C3, 15N) and m+5 glutamine (13C5).
Other tracers can be used by defining their isotope composition (isotopomer) using a space separated list of [<nominal mass>]<element name><element quantity> e.g. the m+3 serine isotopologue with 18O on the hydroxyl and 15N on the amino would be indentified as: [15]N1 [18]O1 (if the element quantity is not supplied it is assumed 1).
Or m+3 serine with uniformly labelled 13C would be indentified as: [13]C3.
These can be specified as `'labels'` in the `params` dictionary, for example:
```
'labels': {
    'ser_m3_1': '[15]N [18]O',
    'ser_m3_2': '[13]C3'
    }
```

These label names must be consistent with those used in the input JSON file that describes the samples.
They are also used to refer back to the peak area ratio: parent/(parent+heavy).
For example using a ratio of 0.4 to 0.6 for ser_m3_1 and 0.5 to 0.7 ser_m3_2:
```
'area_ratio_cutoff': {
    'ser_m3_1': ((0.4, 0.6),),
    'ser_m3_2': ((0.5, 0.7),)
```

Finally, the mass shift for each label must be specified.
This can be input manually like this:
```
'MW_shift': {
    'ser_m3_1': 3.001279886900008,
    'ser_m3_2': 3.010064506020001
```

But can also be calculated in place using the isotope composition:
```
params['MW_shift'] = dict()
for label in params['labels']:
    params['MW_shift'][label] = isotope_obj.isotopes2mass_shift(params['labels'][label])
```


#### A note on natural abundance correction
We do not apply natural abundance correction for the simple reason that we do not know the chemical formula of the unknown RMA traced compounds.
For some RMA traced compounds we of course have their chemical formula and could perform natural abundance correction; however, since the purpose of this method is to identify unknwon metabolites, we chose to treat known and unknown compounds consistently and thus skip natural abundance correction all together.
This choice is not a problem for tracers without substantial bleed over from natural abundance e.g. 13C3, 13C4, {13C3, 15N} etc.
However, tracers with single labels like 15N and 13C, and to a lesser extent tracers with double labels like 15N2, 13C2 and {15N, 13C}, are not as robust.
Prospective users must consider this when designing experiments.



### Input requirements
* Excel file with the output of Compound Discoverer. This file should contain one row per compound and three columns for: "Name", "Molecular Weight", and "RT \[min\]". And then followed by the columns for the peak area per sample e.g.: "Area: KD031720_031720_BD_C01.raw (F1)", "Area: KD031720_031720_BD_C02.raw (F2)", etc. The peak area columns must contain the tag "Area", which it will automatically if the Compound Discoverer output is not altered.
* A JSON file describing each input sample, its type, label etc. See the examples provided in the `projects` folder.


## Dependencies
The following Python packages are required and can be install with `pip` or `conda` commands:
* jupyterlab
* pandas
* HTMLParser
* numpy
* seaborn
* matplotlib

We recommend using an [Anaconda](https://www.anaconda.com/download) Python install which already has many of the above packages installed by default.

