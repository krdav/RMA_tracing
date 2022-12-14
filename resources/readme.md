# Resources for peak processing and filtering


## Python code processing and filtering
Class and function definitions are stored in `helper.py` which works as a container for most Python code.
These are imported and used in the main processing notebook.


## Isotope mass and abundance data from IUPAC

Latest IUPAC information lifted from their "isotopes-matter" webpage:  
https://iupac.org/isotopes-matter/

Isotope masses were downloaded as csv file:  
https://ciaaw.org/data/IUPAC-atomic-masses.csv  
Saved as `IUPAC-atomic-masses.csv`

Abundances were not available as csv, therefore an html table was lifted from:  
https://ciaaw.org/isotopic-abundances.htm  
The table is named "mytable" in the html source and saved as `IUPAC-atomic-abundances.html`, then processed with Pandas `read_html` function.


## Lists of known compounds

A blacklist is stored as `blacklist.tab`

| Description            | MW        | RT    | Polarity | MW_ppm_tol | RT_tol |
|------------------------|-----------|-------|----------|------------|--------|
| GSH in-source fragment | 273.09605 | 9.806 | neg      | 4          | 2      |



A list of compounds known to be labelled is stored as `known_cys_labeled.tab` and `known_gln_labeled.tab`


| Name | Label | Formula                         | RT  |
|------|-------|---------------------------------|-----|
| GSH  | m     | C10 H17 N3 O6 S                 | 9.8 |
| GSH  | m+4   | C7 \[13\]C3 H17 N2 \[15\]N O6 S | 9.8 |



## List of adducts and adduct definitions


List of adducts lifted from the Fiehn lab:  
https://fiehnlab.ucdavis.edu/staff/kind/metabolomics/ms-adduct-calculator/  
Saved as `adducts.tab`  

All the adducts in this table will be searched for in the matching polarity.
To exclude an adduct from the search mark it with the pound sign, like below.

| Adduct_name   |Adduct_mass  | Charge |
|---------------|-------------|--------|
| M+3H          | 1.007276    | 3      |
| M+Na          | 22.989218   | 1      |
| M+H           | 1.007276    | 1      |
| #M+IsoProp+H  | 61.06534    | 1      |




