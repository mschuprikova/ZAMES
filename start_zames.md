<div style="text-align: justify;">

For a run of ZAMES, the following parameters are required: 

* **Control** --- raw results of **control** sequencing by the Sanger method in **.ab1** format. 
This set of files (or one file) represents the results of Sanger sequencing of control 
sample --- sample that does not carry a mutation but have undergone the same experimental 
sample preparation process.

* **Sample** --- raw results of **sample** sequencing by the Sanger method in **.ab1** format. 
This set of files represents the results of Sanger sequencing of the sample under study --- sample containing a mutation at the target site.
* **Reference** --- the nucleotide sequence of a gene, plasmid or genome in **FASTA** format. This file is required for verifying the sequence of chromatograms, as well as for determining whether their orientations correspond.

* **Target Region** --- the nucleotide sequence of the region of interest. This sequence should consist of a short segment of DNA (approximately 50 nucleotides) that contains the target site. During the operation of the ZAMES algorithms, the sequence of the target region enables the alignment of read orientations; it also allows this area to be visualized in the ZAMES execution report, facilitating a comparison between the mutation frequency at the target site and that at non-target positions (i.e., positions lacking modifications).

* **Target Site** --- the position of the target nucleotide in the region of interest (site of expected mutation). This parameter is necessary for determining the position of a target site within a target region and subsequently correlating it with its position as interpreted.

::: {style="border-left: 5px solid #b02725; padding: 10px; background-color: #f5dada; align-items: center; margin-bottom: 10px;"}
<img src="caution_logo.png" alt="icon" width="20" height="25"/> <strong style="color: #b02725">Caution</strong>\
At the moment, ZAMES is undergoing further development. Therefore, to ensure the application functions correctly, a number of requirements must be met:

1) ZAMES only processes files in ABIF format. Therefore, sequencing results in SCF format cannot be processed correctly.

2) The target region sequence must be unique and present in both the reference and the reads, and it must also contain the target site where the mutation is expected to occur. Furthermore, the length of the target region should not be too great (approximately 50 nucleotides is recommended). In addition, the target region sequence must consist exclusively of characters from the deoxyribonucleotide alphabet (A, T, G, C) in uppercase. Otherwise, ZAMES will be unable to process the target region correctly and will terminate its operation.

3) It is important that the numerical value representing the target site does not exceed the length of the target region.

**The developers are actively working to fix bugs. If all instructions are followed precisely, the ZAMES algorithm will process your files correctly. Thank you for your understanding!**
:::

</div>