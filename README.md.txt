# Probabilistic Modeling of the Blood Microbiome in Colorectal Cancer (CRC)

Academic research project developed for the **University of Milan (Università degli Studi di Milano)** in collaboration with **Prof. Marta Rossi** (Department of Clinical Sciences and Community Health - Section of Medical Statistics, Biometry and Epidemiology) and collaborating institutions across Milan.

---

## Project Overview

Colorectal cancer (CRC) compromises the integrity of the intestinal mucosal barrier, promoting the translocation of bacterial DNA into the bloodstream[cite: 1]. High bacterial load and gut permeability correlate with increased risk of CRC.

Analyzing circulating bacterial DNA presents specific statistical challenges:
1. **High Sparsity & High Dimensionality:** Thousands of taxa dominated by rare species.
2. **Uneven Sampling Depth:** Stochastic noise across sequencing runs.
3. **Discrete Count Data:** Abundances constrained to a simplex.

To address these, this project implements a **Dirichlet Multinomial Mixture (DMM)** Bayesian generative model, as proposed by *Holmes et al. (2012)*, to identify microbial dysbiosis signatures in blood samples.

---

## Data Sources & Collaborating Institutions

The clinical data, biological samples, and research context originate from a multi-center collaboration involving the following academic departments and clinical centers in Milan:

1. **Department of Clinical Science and Community Health** (Dipartimento di Eccellenza 2023–2027), University of Milan, Milan, Italy
2. **Department of Biotechnology and Biosciences**, University of Milan-Bicocca, Milan, Italy
3. **Department of Food, Environmental and Nutritional Sciences**, University of Milan, Milan, Italy
4. **Digestive and Interventional Endoscopy Unit**, ASST Grande Ospedale Metropolitano Niguarda, Milan, Italy
5. **Department of Pathophysiology and Transplantation**, University of Milan, Milan, Italy
6. **Gastroenterology and Endoscopy Unit**, Fondazione IRCCS Ca’ Granda Ospedale Maggiore Policlinico, Milan, Italy
7. **Fondazione IRCCS Ca’ Granda Ospedale Maggiore Policlinico**, Milan, Italy

---

## Data Privacy & Confidentiality Notice

The clinical datasets analyzed in this project contain sensitive clinical and genomic information from patients enrolled in the collaborating hospital units listed above. In compliance with strict privacy regulations (GDPR) and institutional ethics requirements, **raw data files are strictly confidential and omitted from this repository**.

---

## Repository Contents

* `Statistica_Medica.R`: R workflow for data processing, Bayesian mixture fitting via Laplace approximation, differential abundance testing, and hybrid clinical models[cite: 1, 2].
* `colon.pdf`: Technical documentation detailing the biological rationale, mathematical formulation of DMMs, and statistical findings.
* `pone.0030126.pdf`: Reference methodology paper (*Holmes et al., 2012*).

---

## Key Findings

* **Optimal Stratification:** Laplace approximation favors $K=2$ partitions, effectively separating healthy controls from CRC-associated dysbiosis states.
* **Anna Karenina Principle (AKP):** Healthy blood microbiome configurations exhibit uniform structure, whereas tumor dysbiosis displays elevated inter-individual variation.
* **Discriminant Taxa:** Key differential variations were identified in *Bacteroides*, *Prevotella*, *Roseburia*, *Faecalibacterium*, and *Ruminococcus*.

---

## Required R Packages

```R
if (!require("haven")) install.packages("haven")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("DirichletMultinomial")) BiocManager::install("DirichletMultinomial")
if (!require("ggplot2")) install.packages("ggplot2")