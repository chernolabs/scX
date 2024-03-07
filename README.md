<h1 align="center">  scX package </h1>

<div align="justify" >
<img align="left" width="25%" src="man/figures/scX.png" title="@irimiodo"> scX is an R package that enables interactive visualization and analysis of single-cell experiment data by creating a Shiny-based application. With scX all aspects of single-cell data can be explored using a large number of different types of plots, such as scatter plots, heatmaps, boxplots, dot and violins plots, etc. All the information associated with cells can be displayed in a customized way: both numerical variables such as logcounts or pseudotime values, and categorical variables such as cell types or sample. One of the main hallmarks of scX is the possibility to plot the main embeddings used for single cell - UMAP, tSNE and PCA - both in 2D and 3D in an interactive way. Thus, embeddings can be rotated, translated and zoomed. But scX is not only a visualization tool, it also allows you to perform different types of analysis on single cell data, such as finding the markers of a cell type or determining the differential genes between two different conditions.
</div>

---

<img src="man/figures/cluster_markers.gif" width="100%" />

<h2>Installation</h2>

<p>scX requires certain packages to be installed. If you do not have these, you can execute the following code:</p>
<div class="sourceCode" id="cb1"><pre class="downlit sourceCode r">
<code class="sourceCode R"><span><span class="fu">BiocManager</span><span class="fu">::</span><span class="fu"><a href="https://rdrr.io/pkg/BiocManager/man/install.html" class="external-link">install</a></span><span class="op">(</span><span class="fu">c</span><span class="op">(</span><span class="st">"SingleCellExperiment"</span>,<span class="st">"scran"</span>,<span class="st">"scater"</span>,<span class="st">"ComplexHeatmap"</span><span class="op">)</span><span class="op">)</span></span></code></pre></div>

<p>scX can be installed from GitHub as follows:</p>
<div class="sourceCode" id="cb1"><pre class="downlit sourceCode r">
<code class="sourceCode R"><span><span class="fu">devtools</span><span class="fu">::</span><span class="fu"><a href="https://remotes.r-lib.org/reference/install_github.html" class="external-link">install_github</a></span><span class="op">(</span><span class="st">"chernolabs/scX"</span><span class="op">)</span></span></code></pre></div>

<p>You can also check the <a href="https://chernolabs.github.io/scX/articles/get_started.html" class="external-link">detailed explanation</a> of what the app can do using the package example data, or interact with our available <a href="http://scx.leloir.org.ar/" class="external-link">app example</a>.</p>
