# Microbiome of Amyotrophic Lateral Sclerosis (ALS) using 16S rRNA

![Stephen Hawking](images/stephen_hawking.png)
<small>*Stephen Hawking*. Photo Credit: [NASA/Paul E. Alers](https://flic.kr/p/6h1t6B).</small>

Amyotrophic lateral sclerosis ([ALS](https://en.wikipedia.org/wiki/Amyotrophic_lateral_sclerosis)) is a neurodegenerative neuromuscular disease that results in the progressive loss of motor neurons that control voluntary muscles. About 20 genes are associated with ALS, most importantly [C9orf72](https://en.wikipedia.org/wiki/C9orf72), which accounts for about 40% of the cases. In addition to genetic risk, environmental factors such as smoking and physical activity represent potential risks. Among the environmental factors is the gut microbiota, which has been shown to contribute and affect mental health, leading to the emerging paradigm of the [gut-brain axis](https://en.wikipedia.org/wiki/Gut%E2%80%93brain_axis) ([GBA](https://en.wikipedia.org/wiki/Gut%E2%80%93brain_axis)). Therefore, [Hertzberg et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33818222/) examined the gut microbiome profiles between ALS patients and their corresponding health caregivers using 16S rRNA.

Here, we are going to reanalyze the 16S dataset generated from [Hertzberg et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33818222/)'s published work. The NCBI BioProject accession number is [PRJNA566436](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA566436), composed of 9 ALS patients and corresponding healthy controls (i.e., a total of 18 samples). The entire workflow is executed within [R](https://www.r-project.org/), mainly using [DADA2](https://github.com/benjjneb/dada2), implemented in the [Bioconductor](https://www.bioconductor.org/) package [`dada2`](https://www.bioconductor.org/packages/release/bioc/html/dada2.html) ([Callahan et al. 2016](https://pubmed.ncbi.nlm.nih.gov/27214047/)). Most of the steps are essentially based on the [DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial.html).

Here is the [analysis](ALS.md) (still work-in-progress) as [GitHub Document](https://rmarkdown.rstudio.com/github_document_format.html), knitted from the [R Markdown](https://rmarkdown.rstudio.com/) source [ALS.Rmd](ALS.Rmd).

**Note**: running on a decent server, the whole analysis takes about 5.5 minutes.
