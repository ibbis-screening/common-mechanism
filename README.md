# An International Common Mechanism for DNA Synthesis Screening

From https://www.frontiersin.org/articles/10.3389/fbioe.2019.00086/full: "As scale drives down cost per base pair, the relatively fixed cost of screening plays a more direct role in overall price. These costs are driven by both customer and sequence screening—commercially-available customer screening solutions still require a great deal of manual review of false positive findings. These false positives create a floor on the possible reduction in labor cost of new customer onboarding. Current sequence screening algorithms are computationally expensive and, given the high false positive rate, the results of sequence screening can be complicated to interpret. These generally require a PhD in bioinformatics both for implementation as well as day to day interpretation of hits. This makes scaling interpretation, in the absence of high-quality sequence annotation, a very expensive proposition." "It is challenging, however, to determine when a custom-built screening system is “good enough”—especially given that the details of each screening implementation remain private to the implementing company. In addition, the recommendations in the Guidance do not specify particular performance metrics in terms of overall sensitivity and specificity or the degree to which sequence alteration or the source of annotation should impact screening results."

## Workflow design
Link to edit: https://docs.google.com/drawings/d/1WTmkvCcyxSCV_KjzTuiacVrUScOpcHy6C59grtiGvTU/edit?usp=sharing
![resources](https://docs.google.com/drawings/d/e/2PACX-1vRQ8uJzbXDgQi68p_S-f6EssH-QgRfuqDhV9QFI4eZRn_CLJJrPbYB8U1n6CWl873G9y-R-q1FdrnNf/pub?w=2570&h=2360)

## User survey
https://docs.google.com/forms/d/1LqzEH3XkFHBtMHNzGW9i4vRU7TyabBnNIDOFQBXa-4s/edit?usp=sharing

## Reading
Search acceleration:
* Data structures based on k-mers for querying large collections of sequencing datasets: https://www.biorxiv.org/content/10.1101/866756v3
* Look into COBS indexing: https://www.biorxiv.org/content/10.1101/2021.03.02.433662v1; https://arxiv.org/abs/1905.09624
  * We have collaborators who could help with this approach
* Look into sketching and containment - mash screen https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1841-x
* Look into ntHash: https://github.com/bcgsc/ntHash

Protein design:
* Discovery, Design, and Structural Characterization of AlkaneProducing Enzymes across the Ferritin-like Superfamily: http://tagkopouloslab.ucdavis.edu/wp-content/uploads/2021/01/discoverydesign.pdf
