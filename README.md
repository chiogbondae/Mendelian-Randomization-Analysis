# Mendelian-Randomization-Analysis
This repository contains the code and results for a genetic study investigating whether Thyroid Stimulating Hormone (TSH) causal effects on stroke are mediated by Atrial Fibrillation (AF). The data and result plots are attached.
source GWAS TSH: Teumer et al. [PMID: 30367059] (TSH: n=54,288)
source GWAS AF: Christophersen et al. [PMID: 28747752] (AF: n=17,931 AF cases + 115,142 controls)
[exclusion.txt](https://github.com/user-attachments/files/24974670/exclusion.txt)


 <img width="3000" height="2400" alt="tsh_stroke_scatter_enhanced" src="https://github.com/user-attachments/assets/10405a83-504e-4982-b74c-dadfb5557d67" />
<img width="3000" height="3600" alt="tsh_stroke_forest_robust" src="https://github.com/user-attachments/assets/7034800f-c99c-46b4-8716-a3533ce58e7c" />
<img width="4200" height="1800" alt="mediation_pathway_diagram" src="https://github.com/user-attachments/assets/e54d9c1f-b07a-4e59-8ff7-27661f744745" />
Results Summary
Pathway	Beta (β)	P-value	Interpretation
TSH → Stroke (Total)	-0.052	3.47×10−2	Significant protective effect
AF → Stroke (Direct)	0.171	5.95×10−30	Major independent risk factor
TSH → Stroke (Direct)	-0.025	0.194	Not significant after adjustment
