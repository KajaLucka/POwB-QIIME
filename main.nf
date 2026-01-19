// Kaja Lucka

nextflow.enable.dsl = 2

params.indir  = "input"       // domyslny katalog z inputem
params.reports  = "reports"   // domyslny katalog z raportami
params.outdir = "output"      // domyslny katalog z outputem


// nazwa mowi sama za siebie - pobieramy dane z interentu i umieszczamy je w katalogu input
process download_data
{
    publishDir params.indir, mode: 'symlink'

    // inputu nie ma, bo dane sa pobierane z internetu

    output:
        path "qiime2-moving-pictures-tutorial"

    script:
    """
        mkdir -p qiime2-moving-pictures-tutorial/emp-single-end-sequences

        wget -O qiime2-moving-pictures-tutorial/sample-metadata.tsv \
            https://data.qiime2.org/2024.10/tutorials/moving-pictures/sample_metadata.tsv

        wget -O qiime2-moving-pictures-tutorial/emp-single-end-sequences/barcodes.fastq.gz \
            https://data.qiime2.org/2024.10/tutorials/moving-pictures/emp-single-end-sequences/barcodes.fastq.gz

        wget -O qiime2-moving-pictures-tutorial/emp-single-end-sequences/sequences.fastq.gz \
            https://data.qiime2.org/2024.10/tutorials/moving-pictures/emp-single-end-sequences/sequences.fastq.gz
    """
}


// importujemy pobrane sekwencje do qiime
process import_sequences_to_qiime
{
    publishDir params.outdir, mode: 'symlink'

    conda "${projectDir}/envs/powb-qiime.yml"  // skorzystaj ze srodowiska z qiime

    input:
        path seseq_dir  // folder z plikami fastq.gz

    output:
        path "emp-single-end-sequences.qza"

    script:
    """
        qiime tools import \
            --type EMPSingleEndSequences \
            --input-path ${seseq_dir} \
            --output-path emp-single-end-sequences.qza
    """
}


process demultiplexing
{
    publishDir params.outdir, mode: 'symlink'

    conda "${projectDir}/envs/powb-qiime.yml"  // skorzystaj ze srodowiska z qiime

    input:
        path ese_seq
        path metadata

    output:
        tuple path("demux.qza"), path("demux-details.qza")

    script:
    """
        qiime demux emp-single \
            --i-seqs ${ese_seq} \
            --m-barcodes-file ${metadata} \
            --m-barcodes-column barcode-sequence \
            --o-per-sample-sequences demux.qza \
            --o-error-correction-details demux-details.qza
    """
}


process demultiplexing_report
{
    publishDir params.reports, mode: 'symlink'

    conda "${projectDir}/envs/powb-qiime.yml"  // skorzystaj ze srodowiska z qiime

    input:
        path demux

    output:
        path "demux.qzv"

    script:
    """
        qiime demux summarize \
            --i-data ${demux} \
            --o-visualization demux.qzv
    """
}

// QC - wybor jednej z dwoch opcji
process qc_dada2
{
    publishDir params.outdir, mode: 'symlink'

    conda "${projectDir}/envs/powb-qiime.yml"  // skorzystaj ze srodowiska z qiime

    input:
        path demux

    output:
        tuple path("rep-seqs.qza"), path("table.qza"), path("stats-dada2.qza")

    script:
    """
        qiime dada2 denoise-single \
            --i-demultiplexed-seqs ${demux} \
            --p-trim-left 0 \
            --p-trunc-len 120 \
            --o-representative-sequences rep-seqs.qza \
            --o-table table.qza \
            --o-denoising-stats stats-dada2.qza
    """
}

process qc_dada2_vis
{
    publishDir params.reports, mode: 'symlink'

    conda "${projectDir}/envs/powb-qiime.yml"  // skorzystaj ze srodowiska z qiime

    input:
        path qc_stats

    output:
        path "stats-dada2.qzv"

    script:
    """
        qiime metadata tabulate \
            --m-input-file ${qc_stats} \
            --o-visualization stats-dada2.qzv
    """
}


process qc_deblur_step1
{
    publishDir params.outdir, mode: 'symlink'

    conda "${projectDir}/envs/powb-qiime.yml"  // skorzystaj ze srodowiska z qiime

    input:
        path demux

    output:
        tuple path("demux-filtered.qza"), path("demux-filter-stats.qza")

    script:
    """
        qiime quality-filter q-score \
            --i-demux ${demux} \
            --o-filtered-sequences demux-filtered.qza \
            --o-filter-stats demux-filter-stats.qza
    """
}

process qc_deblur_step2
{
    publishDir params.outdir, mode: 'symlink'

    conda "${projectDir}/envs/powb-qiime.yml"  // skorzystaj ze srodowiska z qiime

    input:
        path filtered_demux

    output:
    tuple path("rep-seqs.qza"), path("table.qza"), path("deblur-stats.qza")

    script:
    """
        qiime deblur denoise-16S \
            --i-demultiplexed-seqs ${filtered_demux} \
            --p-trim-length 120 \
            --o-representative-sequences rep-seqs.qza \
            --o-table table.qza \
            --p-sample-stats \
            --o-stats deblur-stats.qza
    """
}

process qc_deblur_vis
{
    publishDir params.reports, mode: 'symlink'

    conda "${projectDir}/envs/powb-qiime.yml"  // skorzystaj ze srodowiska z qiime

    input:
        path filtered_stats
        path stats

    output:
        tuple path("demux-filter-stats.qzv"), path("deblur-stats.qzv")

    script:
    """
        qiime metadata tabulate \
            --m-input-file ${filtered_stats} \
            --o-visualization demux-filter-stats.qzv

        qiime deblur visualize-stats \
            --i-deblur-stats ${stats} \
            --o-visualization deblur-stats.qzv
    """
}


workflow
{
    download_data_ch = download_data()

    input_subdir_ch = download_data_ch.map {it + "/emp-single-end-sequences"}
    import_sequences_to_qiime_ch = import_sequences_to_qiime(input_subdir_ch)  // out: emp-single-end-sequences.qza

    metadata_ch = download_data_ch.map {it + "/sample-metadata.tsv"}
    demultiplexing_ch = demultiplexing(import_sequences_to_qiime_ch, metadata_ch)

    demux_ch = demultiplexing_ch.map {it[0]}
    demultiplexing_report_ch = demultiplexing_report(demux_ch)

    if (params.qc == "dada2")
    {
        qc_ch = qc_dada2(demux_ch)

        qc_stats_ch = qc_ch.map{it[2]}
        qc_vis_ch = qc_dada2_vis(qc_stats_ch)
    }
    else if (params.qc == "deblur")
    {
        qc1_ch = qc_deblur_step1(demux_ch)

        qc_filtered_ch = qc1_ch.map{it[0]}
        qc2_ch = qc_deblur_step2(qc_filtered_ch)

        qc_filtered_stats_ch = qc1_ch.map{it[1]}
        qc_stats_ch = qc2_ch.map{it[2]}
        qc_vis_ch = qc_deblur_vis(qc_filtered_stats_ch, qc_stats_ch)
    }
    else
    {
        error "Nieprawidlowa opcja qc = ${params.qc} (dozwolone: dada2, deblur)"
    }
}
