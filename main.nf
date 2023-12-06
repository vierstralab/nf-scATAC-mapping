#!/usr/bin/env nextflow

params.conda = "$moduleDir/environment.yml"

process find_unique_barcodes {

    publishDir "${params.outdir}/barcodes_maps"
    tag "${file_id}"

    input:
        tuple val(file_id), path(fragment_file)
    
    output:
        tuple val(file_id), path(name)

    script:
    name = "${file_id}.unique_barcodes.map"
    """
    cut -f4 ${fragment_file} | sort | uniq > ${name}
    """
}


process split_masterlist_in_chunks {

    publishDir "${params.outdir}/index_chunks"

    output:
        path "${prefix}*"

    script:
    prefix = "chunk_"
    """
    split -l ${params.chunk_size} \
        --additional-suffix=.bed \
        ${params.index_file} \
        ${prefix} 
    """
}


process intersect_with_chunk {
    conda params.conda
    tag "${file_id}:${chunk_id}"
    //scratch true

    input:
        tuple val(file_id), path(barcodes_map), val(chunk_id), path(index_chunk), path(fragment_file)
    
    output:
        tuple val(chunk_id), path(barcodes_map), path(name)

    script:
    name = "${file_id}.${chunk_id}.barcodes.npz"
    """
    bedtools intersect \
        -a ${index_chunk} \
        -b ${fragment_file} \
        -wa -wb \
        | cut -f4,16 > tmp.txt
    
    cut -f4 ${index_chunk} > index_mapping.txt

    python3 $moduleDir/bin/convert_to_sparse_matrix.py \
        tmp.txt \
        ${barcodes_map} \
        index_mapping.txt \
        ${name}
    """
}


process merge_chunks_horizontally {
    conda params.conda
    publishDir "${params.outdir}/matrix_chunks"

    input:
        path sparse_matrices
    
    output:
        path name
    
    script:
    name = ""
    """
    python3 $moduleDir/bin/merge_chunks.py \
        ${name} \
        ${sparse_matrices}
    """

}


workflow map2Index {
    take:
        fragment_files
    main:
        chunks = split_masterlist_in_chunks() 
            | flatten()
            | map(it -> tuple(it.simpleName, it)) // chunk_id, chunk

        out = fragment_files // id, fragment_file
            | find_unique_barcodes // id, barcodes_map
            | combine(chunks) // id, barcodes_map, chunk_id, chunk
            | combine(fragment_files, by:0) // id, barcodes_map, chunk_id, chunk, fragment_file
            | intersect_with_chunk 
            //| groupTuple()
            //| merge_chunks_horizontally()
    emit:
        out
}


workflow {
    Channel.fromPath(params.samples_file)
        | splitCsv(header: true, sep: '\t')
        | map(row -> tuple(row.fragment_file_id, file(row.fragment_file)))
        | map2Index
}
