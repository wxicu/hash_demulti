#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { preprocessing_hashing } from './modules/preprocess'
include { multiseq_hashing } from './modules/multiseq'
include { htodemux_hashing } from './modules/htodemux'
include { hash_solo_hashing } from './modules/hashsolo'
include { hashedDrops_hashing } from './modules/hashedDrops'
include { demuxem_hashing } from './modules/demuxem'
include { solo_hashing } from './modules/solo'

process summary{
    publishDir "$params.outdir/compare", mode: 'copy'
    input:
        val demuxem_result
        val hashsolo_result
        val htodemux_result
        val multiseq_result
        val hashedDrops_result
    
    output:
        path '*.csv'

    script:
        def demuxem_files = ""
        def htodemux_files = ""
        def hashsolo_files = ""
        def multiseq_files = ""
        def hashedDrops_files = ""
        
        if (demuxem_result != "no_result"){
            demuxem_files = "--demuxem "
            for(r : demuxem_result) {
                demuxem_files  = demuxem_files + r + ":"
            }
        }
        if (hashsolo_result != "no_result"){
            hashsolo_files = "--hashsolo "
            for(r : hashsolo_result) {
                hashsolo_files = hashsolo_files + r + ":"
            }
        }
        if (htodemux_result != "no_result"){
            htodemux_files = "--htodemux "
            for(r : htodemux_result) {
                htodemux_files = htodemux_files + r + ":"
            }
        }
        if (multiseq_result != "no_result"){
            multiseq_files = "--multiseq "
            for(r : multiseq_result) {
                multiseq_files = multiseq_files + r + ":"
            }
        }
        if (hashedDrops_result != "no_result"){
            hashedDrops_files = "--hashedDrops "
            for(r : hashedDrops_result) {
                hashedDrops_files = hashedDrops_files + r + ":"
            }
        }
        
        """
        summary.R $demuxem_files $htodemux_files $multiseq_files $hashedDrops_files $hashsolo_files
        """
}



workflow{

    if ((params.htodemux == "True" & params.htodemux_preprocess != "False")| \
       (params.multiseq == "True" & params.multiseq_preprocess != 'False')){
        preprocessing_hashing()
    }
    
    if (params.htodemux == "True"){
        rdsobj = params.htodemux_preprocess == 'True'? preprocessing_hashing.out: (params.htodemux_preprocess == 'False'? Channel.from(params.rdsObj_htodemux) : preprocessing_hashing.out.mix(Channel.from(params.rdsObj_htodemux)))
        htodemux_hashing(rdsobj)
        htodemux_out = htodemux_hashing.out
    }
    else{
        htodemux_out = channel.value("no_result")
    }
    
    if (params.multiseq == "True"){
        rdsobj = params.multiseq_preprocess == 'True'? preprocessing_hashing.out: (params.multiseq_preprocess == 'False'? Channel.from(params.rdsObj_multi) : preprocessing_hashing.out.mix(Channel.from(params.rdsObj_multi)))
        multiseq_hashing(rdsobj)
        multiseq_out = multiseq_hashing.out
    }
    else{
        multiseq_out = channel.value("no_result")
    }
    
    if (params.hashsolo == "True"){
        hash_solo_hashing()
        hashsolo_out = hash_solo_hashing.out
    }
    else{
        hashsolo_out = channel.value("no_result")
    }
    
    if (params.demuxem == "True"){
        demuxem_hashing()
        demuxem_out = demuxem_hashing.out
    }
    else{
        demuxem_out = channel.value("no_result")
    }
    
    if (params.hashedDrops == "True"){
        hashedDrops_hashing()
        hashedDrops_out = hashedDrops_hashing.out
    }
    else{
        hashedDrops_out = channel.value("no_result")
    }
    
    if (params.solo == "True"){
        solo_hashing()
        solo_out = solo_hashing.out
    }
    else{
        solo_out = channel.value("no_result")
    }
    
    summary(demuxem_out, hashsolo_out, htodemux_out, multiseq_out, hashedDrops_out)

}
