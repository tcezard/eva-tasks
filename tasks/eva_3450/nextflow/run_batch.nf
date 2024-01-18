#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
    Process a chunk of merge or split candidates

    Inputs:
            --source_duplicate_file     File containing list of duplicate (split or merge)
            --clustering_props          Properties file for the Merge resolution (Clustering CLUSTER_UNCLUSTERED_VARIANTS_JOB)
            --assembly_accession        Target assembly where the merge and split should be detected and corrected
            --instance_id               Instance id to run clustering
            --output_dir                Directory where the output log will be copied 
    """
}

params.source_duplicate_file = null
params.clustering_props = null
params.assembly_accession = null
params.instance_id = null
params.output_dir = null
params.write_properties_file=''
params.read_properties_file=''
params.load_merge_groovy=''
params.load_split_groovy=''
params.clustering_jar=''

if (!params.source_duplicate_file || !params.clustering_props || !params.instance_id || !params.output_dir || !params.assembly_accession){
    if (!params.source_duplicate_file) log.warn('Provide a file containing the list of duplicates to resolve using --source_duplicate_file')
    if (!params.clustering_props) log.warn('Provide a file with the properties for the Clustering job using --clustering_props')
    if (!params.assembly_accession) log.warn('Provide an assembly accession where the candidate have been detected using --assembly_accession')
    if (!params.instance_id) log.warn('Provide an instance id using --instance_id')
    if (!params.output_dir) log.warn('Provide an output directory using --output_dir')
    exit 1, helpMessage()
}

workflow load_split {
    load_split_candidate_variants(params.source_duplicate_file)
    remediate_split_and_merge(load_split_candidate_variants.out.candidate_load_log)
}

workflow load_merge {
    load_merge_candidate_variants(params.source_duplicate_file)
    remediate_split_and_merge(load_merge_candidate_variants.out.candidate_load_log)
}


process load_merge_candidate_variants {
    memory '8 GB'

    input:
    path source_duplicate_file

    output:
    path "${source_duplicate_file}_merged_candidate.log", emit: candidate_load_log

    publishDir "$params.output_dir", overwrite: true, mode: "copy"

    script:
    groovy_project_dir=file(params.load_merge_groovy).getParent().getParent().getParent().getParent().getParent()
    """
    run_groovy_script.sh $groovy_project_dir $params.load_merge_groovy -assemblyAccession=$params.assembly_accession -devPropertiesFile=$params.write_properties_file  -prodPropertiesFile=$params.read_properties_file -rsMergeCandidateFile=${source_duplicate_file.toRealPath()} > ${source_duplicate_file}_merged_candidate.log 2>&1
    """
}

process load_split_candidate_variants {
    memory '8 GB'

    input:
    path source_duplicate_file

    output:
    path "${source_duplicate_file}_split_candidate.log", emit: candidate_load_log

    publishDir "$params.output_dir", overwrite: true, mode: "copy"

    script:
    groovy_project_dir=file(params.load_split_groovy).getParent().getParent().getParent().getParent().getParent()
    """
    run_groovy_script.sh $groovy_project_dir $params.load_split_groovy -assemblyAccession=$params.assembly_accession -devPropertiesFile=$params.write_properties_file  -prodPropertiesFile=$params.read_properties_file -rsSplitCandidateFile=${source_duplicate_file.toRealPath()} > ${source_duplicate_file}_split_candidate.log 2>&1
    """
}


process remediate_split_and_merge {
    memory '8 GB'
    clusterOptions "-g /accession/instance-${params.instance_id}"

    input:
    path candidate_load_log

    output:
    path "${candidate_load_log.getBaseName()}_clustering.log", emit: clustering_log_filename
    path "${candidate_load_log.getBaseName()}_rs_report.txt", optional: true, emit: rs_report_filename

    publishDir "$params.output_dir", overwrite: true, mode: "copy"

    """
    java -Xmx8G -jar $params.clustering_jar \
        --spring.config.location=file:${params.clustering_props} \
        --spring.batch.job.names=CLUSTER_UNCLUSTERED_VARIANTS_JOB \
        --parameters.rsReportPath=${candidate_load_log.getBaseName()}_rs_report.txt \
        --parameters.assemblyAccession=$params.assembly_accession \
        --accessioning.instanceId=instance-${params.instance_id} \
        --parameters.remappedFrom=Dummy \
        --parameters.projects=Dummy \
        --parameters.projectAccession=Dummy \
        --parameters.vcf=Dummy \
        > ${candidate_load_log.getBaseName()}_clustering.log
    """
}
