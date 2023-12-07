#!/usr/bin/env nextflow

params.conda = "$moduleDir/environment.yml"

process find_unique_barcodes {

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

process create_dhs_map {

    publishDir params.outdir
    output:
        path name


    script:
    name = "index_mapping.txt"
    """
    cut -f4 ${params.index_file} > ${name}
    """
}


process intersect_with_index {
    conda params.conda
    tag "${file_id}"
    label "med_mem"
    scratch true

    input:
        tuple val(file_id), path(barcodes_map), path(dhs_map), path(fragment_file)
    
    output:
        tuple path(name), path(barcodes_map)

    script:
    name = "${file_id}.barcodes.npz"
    """
    sort-bed ${fragment_file} | bedtools intersect \
        -a ${params.index_file} \
        -b stdin \
        -wa -wb -sorted \
        | cut -f4,16 > tmp.txt
    

    python3 $moduleDir/bin/convert_to_sparse_matrix.py \
        tmp.txt \
        ${barcodes_map} \
        ${dhs_map} \
        ${name}
    """
}


process merge_chunks {
    conda params.conda
    publishDir params.outdir
    label "high_mem"

    input:
        path matrices_and_maps
    
    output:
        path "${prefix}*"
    
    script:
    prefix = "all_barcodes"
    """
    python3 $moduleDir/bin/merge_chunks.py \
        ${params.samples_file} \
        ${prefix}
    """

}


workflow map2Index {
    take:
        fragment_files
    main:
        out = fragment_files // id, fragment_file
            | find_unique_barcodes // id, barcodes_map
            | combine(
                create_dhs_map()
            ) // id, barcodes_map, dhs_map
            | combine(fragment_files, by:0) // id, barcodes_map, dhs_map, fragment_file
            | intersect_with_index // barcodes_map, sparse_matrix
            | collect(flat: true, sort: true)
            | merge_chunks // matrices_and_maps
    emit:
        out
}


workflow {
    Channel.fromPath(params.samples_file)
        | splitCsv(header: true, sep: '\t')
        | map(row -> tuple(row.fragment_file_id, file(row.fragment_file)))
        | map2Index
}
