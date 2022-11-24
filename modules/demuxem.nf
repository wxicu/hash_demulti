process DEMUXEM{
    publishDir params.outdir/demuxem, mode:'copy'
    label "demuxem"
    input:
        path rna_raw 
        path hto_mat_folder
	val gz
	val n_threads
	val genome
        val alpha
        val min_num_genes
        val min_num_umis
        val min_signal_hashtag
        val random_state
	val generate_diagnostic_plots
	val generate_gender_plot 
        val output_demux
    output:
	path 'matrix_filtered.csv'
        path '*_demux.zarr'
	path '*.out.demuxEM.zarr'
	path '*.ambient_hashtag.hist.pdf'
	path '*.background_probabilities.bar.pdf'
	path '*.real_content.hist.pdf'
	path '*.rna_demux.hist.pdf'
	path '*.gene_name.violin.pdf'
        
    script:
	def generateDiagnosticPlots =  generate_diagnostic_plots == 'True' ? " -generate-diagnostic-plots" : ' '
	def generateGenderPlot =  generate_gender_plot != 'False' ? " -generate-gender-plot ${generate_gender_plot}" : ' '
        """
        python mat2csv.py --hto_mat_folder $hto_mat_folder --gz $gz
	demuxEM -p ${n_threads} -alpha-on-samples ${alpha} -min-num-genes ${min_num_genes} -min-num-umis ${min_num_umis} \
		-min-signal-hashtag ${min_signal_hashtag} -random-state ${random_state} ${generateDiagnosticPlots} \
		${generateGenderPlot} ${rna_raw} matrix_filtered.csv ${output_demux}

        """

}




workflow demuxem_hashing{
  main:
	rna_raw = Channel.fromPath(params.rna_raw)
        hto_mat_folder = Channel.fromPath(params.hto)
  	gz = Channel.from(params.gz)
        n_threads = Channel.from(params.n_threads)
	genome = Channel.from(params.genome)
	alpha = Channel.from(params.alpha)
        min_num_genes = Channel.from(params.min_num_genes)
        min_num_umis = Channel.from(params.min_num_umis)
        min_signal_hashtag = Channel.from(params.min_signal_hashtag)
        random_state = Channel.from(params.random_state)
        generate_diagnostic_plots = Channel.from(params.generate_diagnostic_plots)
        generate_gender_plot = Channel.from(params.generate_gender_plot)
        output_demux = Channel.from(params.output_demux)

	DEMUXEM(rna_raw, hto_mat_folder, gz, n_threads, genome, alpha, min_num_genes, min_num_umis, min_signal_hashtag, random_state, generate_diagnostic_plots, generate_gender_plot, output_demux)
  
  emit:
	DEMUXEM.out


}
