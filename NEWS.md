# SpatialEcoTyper 0.0.7

* First release of SpatialEcoTyper.
  * SpatialEcoTyper: SE discovery from a single sample
  * MultiSpatialEcoTyper: Integrative analysis of SEs from multiple samples
  * IntegrateSpatialEcoTyper: Integrative analysis of SEs from multiple samples
  * RecoverSE: SE recovery
  * nmfClustering: NMF clustering
  * SpatialView: draw spatial map of the tissue
  * HeatmapView: draw a heatmap
  * CreatePseudobulks: create pseudobulk mixtures
  * NMFGenerateW: train an NMF model for SE deconvolution from bulk expression profiles
  * NMFGenerateWList: train cell-type specific NMF model for SE recovery

# SpatialEcoTyper 0.0.6
* NMFpredict: warning when there are limited genes overlap with the model.

# SpatialEcoTyper 0.0.5
  * Reduce memory usage by computing distance of spatial neighbors and not all pairwise distances

# SpatialEcoTyper 0.0.4
  * Reduce memory usage by replacing Reduce() with loops
  * Add seeds to nmfClustering
  * Re-organize the documentation for integrative analysis
  * Test and refine all documentations
  
# SpatialEcoTyper 0.0.3
  * Test Seurat v4.2, v4.4 and v5 for the analysis and add related notes: they lead to different embedding and clustering results, but show high consistency (ARI=0.7) for the demo.
  * Add more figures to the output directory for integrative analysis
  * Add hints about the training of SE recovery model: the demo is less robust due to limited number of cells used. The training data should be as comprehensive as possible.
  * Add minibatch option to SNF2 function to reduce memory usage
  * Add hints about memory usage and parallel processing

# SpatialEcoTyper 0.0.2
  * Add `filter.region.by.celltypes` option to SpatialEcoTyper and MultiSpatialEcoTyper

# SpatialEcoTyper 0.0.1
  * Second release
  * Update SE recovery from Visium data
  * Add option `dropcell` to SpatialEcoTyper function to drop NA regions

