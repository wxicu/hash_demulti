process PREPROCESS{
    publishDir params.outdir/preprocess, mode:'copy'
    label "seurat_process"
    input:
        val rdsObject
        val umi_counts
        val hto_matrix
        val ndelim
        val selection_method
        val number_features
        val assay
        val margin
        val normalisation_method
        val preprocess_out
    output:
        file '*.rds'

    script:
    """
        Rscript pre_processing.R  --rdsObject $rdsObject --fileUmi $umi_counts --fileHto $hto_matrix --ndelim $ndelim \
                                  --selectMethod $selection_method --numberFeatures $number_features --assay $assay \
                                  --margin $margin --normalisationMethod $normalisation_method --OutputFile $preprocess_out
    """


}


workflow preprocessing{
    main:
        rdsObject = Channel.from(params.rdsObject)
        umi_matrix = Channel.fromPath(params.umi_count)
        hto_matrix =  Channel.fromPath(params.hto_mat)
        sel_method = Channel.from(params.selection_method)
        ndelim = Channel.from(params.ndelim)
        n_features = Channel.from(params.number_features)
        assay = Channel.from(params.assay)
        margin = Channel.from(params.margin)
        norm_method = Channel.from(params.normalisation_method)
        out_file = Channel.from(params.preprocessOut)
        PREPROCESS(rdsObject, umi_matrix, hto_matrix, ndelim, sel_method, n_features, assay, margin, norm_method, out_file)

}
  