process MULTI_SEQ{
    publishDir params.outdir/multiseq, mode:'copy'
    
    input:
    val rdsObject
    val quantile_multi
    val autoThresh
    val maxiter
    val qrangeFrom
    val qrangeTo
    val qrangeBy
    val verbose
    val assay
    val objectOutMulti
    val assignmentOutMulti


    output:
        path '*.rds'
        path '*_assignment_multiseq_demux.csv'

    script:
    """
    
    Rscript MultiSeq.R --rdsObject $rdsObject --assay $assay --quantile_multi $quantile_multi --autoThresh $autoThresh --maxiter $maxiter --qrangeFrom $qrangeFrom --qrangeTo $qrangeTo --qrangeBy $qrangeBy --verbose $verbose  --objectOutMulti $objectOutMulti --assignmentOutMulti $assignmentOutMulti
    """
}



workflow multiseq_hashing{
   main:
        rdsObject =  Channel.fromPath(params.rds_object)
  	quantile_multi = Channel.from(params.quantile_multi)
  	autoThresh = Channel.from(params.autoThresh)
  	maxIter = Channel.from(params.maxiter)
  	qrangeFrom = Channel.from(params.qrangeFrom)
  	qrangeTo = Channel.from(params.qrangeTo)
  	qrangeBy = Channel.from(params.qrangeBy)
  	verbose = Channel.from(params.verbose)
	assay = Channel.from(params.assay)
  	out_multi = Channel.from(params.objectOutMulti)
  	assignment_multi = Channel.from(params.assignmentOutMulti)
        MULTI_SEQ(rdsObject, quantile_multi, autoThresh, maxIter, qrangeFrom, qrangeTo, qrangeBy, verbose, assay, out_multi, assignment_mutli)
   emit:
	MULTI_SEQ.out


} 