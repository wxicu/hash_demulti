nextflow.enable.dsl=2
process HTODEMUX{
    publishDir "$params.outdir/htodemux", mode: 'copy'
    label "seurat_process"
    input:
        path seurat_object
        val quantile_hto
        val kfunc
        val nstarts
        val nsamples
        val seed
        val init
        val assignmentOutHTO
        val objectOutHTO

    output:
        path '*.rds'
        path '*_classification_htodemux.csv'
        path '*_assignment_htodemux.csv'

    script:

        """
            Rscript HTODemux.R --seuratObjectPath $seurat_object --quantile $quantile_hto --kfunc $kfunc --nstarts $nstarts --nsamples $nsamples --seed $seed --init $init --assignmentOutHTO $assignmentOutHTO --objectOutHTO $objectOutHTO
        """


}

process HTO_VISUALISATION{

    publishDir path: "$projectDir/htodemux", mode:'copy'
    label "seurat_process"
    
    input:
    path result_object
    val assay
    //Ridge plot params
    val ridgePlot
    val ridgeNCol
    //Scatter features params
    val featureScatter
    val scatterFeat1
    val scatterFeat2
    //Violin plot params
    val vlnplot
    val vlnFeatures
    val vlnLog
    //tSNE
    val tsne
    val tseIdents
    val tsneInvert
    val tsneVerbose
    val tsneApprox
    val tsneDimMax
    val tsePerplexity
    //Heatmap
    val heatmap
    val heatmapNcells

    output:
    file '*.png'

    script:
    """
        Rscript HTODemux-visualisation.R --HashtagPath $result_object --assay $assay --ridgePlot $ridgePlot --ridgeNCol $ridgeNCol --featureScatter $featureScatter --scatterFeat1 $scatterFeat1 --scatterFeat2 $scatterFeat2 --vlnplot $vlnplot --vlnFeatures $vlnFeatures --vlnLog $vlnLog --tsne $tsne --tseIdents $tseIdents --tsneInvert $tsneInvert --tsneVerbose $tsneVerbose --tsneApprox $tsneApprox --tsneDimMax $tsneDimMax --tsePerplexity $tsePerplexity --heatmap $heatmap --heatmapNcells $heatmapNcells
    """


}


workflow HTODemux_hashing{
    take:
	seurat_object
    main:
        seurat_object
        quantile_hto = Channel.from(params.quantile_hto)
  	kfunc = Channel.from(params.kfunc)
  	n_starts = Channel.from(params.nstarts)
  	n_samples = Channel.from(params.nsamples)
  	seed = Channel.from(params.seed)
  	init = Channel.from(params.init)
  	object_hto = Channel.from(params.objectOutHTO)
  	assignment_hto = Channel.from(params.assignmentOutHTO)



        HTODEMUX(seurat_object, quantile_hto, kfunc, nstarts, nsamples, seed, init, assignmentOutHTO, objectOutHTO)
        if visualize_HTODEMUX == TRUE{
           HTO_VISUALISATION(HTODEMUX.out[0])
           
}
    emit:
        HTODEMUX.out


}