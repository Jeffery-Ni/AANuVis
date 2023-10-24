AANuVis
==========================================================================================================================
***A semi-automatic amino acid usage and di- and tetra-nucleotide usage frequency calculator and UMAP visualizer***

# Step one
Use # AANuVis-count_local.py # for local fasta file AAU and di- and tetra-nucleotide usage count and convert into percentaged format

   Usage: ```python AANuVis-count_local.py fasta_file.fasta``` or ```nohup AANuVis-count_local.py fasta_file.fasta &``` for none-verbose
   
   **Dependency: seqkit, prodigal, python**

# Step two
Then use # AANuVis-Visualize.py # for visualization. Parameters are integrated in the code, feel free to adjust for individualized and customized uses. :-) 

  Usage: You have to change ```csv_files =```  ```colors =``` and ```dot_alphas =``` in the code. For visualization customization like ```grid_interval``` ```grid_alpha```, ```line_width```,```x_tick_fontsize```,```y_tick_fontsize```, and much much more, see annotations in the code.  
  
  Calculation Parameters are in umap.UMAP, this is where you need to know how UMAP works, see 10.1038/nbt.4314 for introduction, or just any youtuber's video

  **Dependency: numpy, matplotlib.pyplot, umap-learn**
