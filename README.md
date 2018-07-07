# clusteringAbstracts

`Shiny` app to cluster abstracts of scienctific journals.

A 2 cols (id and abstract in this order) csv file, or a [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/) retrieved XML file are applicable.

## Method
Distance is cosine distance.
Hierarchial clustering is calculated by Ward method.

## Output
- **Dendrogram**: showing cluster number and key word branch
- **Top word table**: Most 10 specific words of each cluster
- **Key word(s) density table**: Cluster size in percentage, and key word(s) density in each cluster

## Usage
1. Load data file from "Select file:".
(2. Optional, but highly recommended: load stop-words file from "Stop Words File:".)
3. Click "Cluster" to plot dendrogram.
4. Decide cluster number, and click "Cluster" again to draw dendrogram with cluster rectangle.
5. Type key words into "Keywords:" area.
6. Click "Label" to make key word density table.
