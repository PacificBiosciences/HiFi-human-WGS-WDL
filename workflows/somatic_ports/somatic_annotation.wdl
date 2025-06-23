version 1.0



## task definitons 
task vep_annotate {
    input {
        File input_vcf
        File? vep_cache
        File ref_fasta
        File ref_fasta_index
        Int threads
    }

    Float file_size = ceil(size(input_vcf, "GB") + size(vep_cache, "GB") + size(ref_fasta, "GB") + size(ref_fasta_index, "GB"))

    command <<<
        set -euxo pipefail

        mkdir -p vep_data/
        # If vep_cache is not provided, fail
        if [ ! -f ~{vep_cache} ]; then
            echo "VEP cache file not found. Please provide a valid cache file."
            exit 1
        fi

        vep --help

        tar -xzvf ~{vep_cache} -C vep_data/
        vep \
            --cache \
            --offline \
            --dir vep_data/ \
            --fasta ~{ref_fasta} \
            --format vcf \
            --fork ~{threads} \
            --species homo_sapiens \
            --assembly GRCh38 \
            --symbol \
            --hgvs \
            --refseq \
            --check_existing \
            --vcf \
            --pick \
            --flag_pick_allele_gene \
            --everything \
            --compress_output bgzip \
            -i ~{input_vcf} \
            -o ~{sub(basename(input_vcf), "\\.vcf.gz$", "")}.vep.vcf.gz

        # Delete cache after annotation
        rm -rf vep_data/

    >>>

    output {
        File vep_annotated_vcf = sub(basename(input_vcf), "\\.vcf.gz$", "") + ".vep.vcf.gz"
    }

    runtime {
        docker: "ensemblorg/ensembl-vep@sha256:e7612ab7c2923f2b9a78592b939e74874cd29f7494d70ee7135c8303841b03a8"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task annotsv {
    input {
        File sv_vcf
        File sv_vcf_index
        File? annotsv_cache
        Int threads
    }

    Float file_size = ceil(size(sv_vcf, "GB") + size(annotsv_cache, "GB"))

    command <<<
        set -euxo pipefail

        # Process VCF to move BND alt to INFO field
        # Check if it ends in gz. Unzip it if it's the case
        if [[ ~{sv_vcf} == *.gz ]]; then
            gunzip -c ~{sv_vcf} > tmp.vcf
        else
            cp ~{sv_vcf} tmp.vcf
        fi
        
        awk -F'\t' -v OFS='\t' '
        { if (NR==2) {
                print "##INFO=<ID=SV_ALT,Number=1,Type=String,Description=\"Square bracketed notation for BND event\">"
            }
        }
        {
            if ($0 ~ /^#/) {
                print $0;  # Print header lines as is
            } else {
                if ($8 ~ /SVTYPE=BND/) {
                    $8 = $8 ";SV_ALT=" $5;  # Append ALT to INFO as SV_ALT
                    $5 = "<BND>";  # Change ALT to <BND>
                }
                print $0;
            }
        }' tmp.vcf > tmp_processed.vcf

        mkdir -p annotsv_cache_dir
        # If annotsv_cache is not provided, fail
        if [ ! -f ~{annotsv_cache} ]; then
            echo "AnnotSV cache file not found. Please provide a valid cache file."
            exit 1
        fi

        AnnotSV --version

        tar -xzvf ~{annotsv_cache} -C annotsv_cache_dir/

        AnnotSV \
            -annotationsDir annotsv_cache_dir/AnnotSV/ \
            -SVinputFile tmp_processed.vcf \
            -outputDir . \
            -outputFile ~{sub(basename(sv_vcf), "\\.vcf.gz$", "")}.annotsv.tsv \
            -SVinputInfo 1 \
            -genomeBuild GRCh38

        # Delete cache after annotation
        rm -rf annotsv_cache_dir/
        (head -1 ~{sub(basename(sv_vcf), "\\.vcf.gz$", "")}.annotsv.tsv  | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$5,$6,$14,$26,$54,$81,$104,$106,$107,$118}' && 
         tail -n +2 ~{sub(basename(sv_vcf), "\\.vcf.gz$", "")}.annotsv.tsv | 
         awk -F '\t' -v OFS='\t' '
         {
            split($107, classifications, /;/);
            max_priority = 0;
            for (i in classifications) {
                gsub(/^[ \t]+|[ \t]+$/, "", classifications[i]);
                if (classifications[i] == "Definitive") curr_priority = 4;
                else if (classifications[i] == "Strong") curr_priority = 3;
                else if (classifications[i] == "Moderate") curr_priority = 2;
                else if (classifications[i] == "Supportive") curr_priority = 1;
                else curr_priority = 0;
                if (curr_priority > max_priority) max_priority = curr_priority;
            }
            print max_priority, ($118+0), $1,$2,$3,$5,$6,$14,$26,$54,$81,$104,$106,$107,$118
         }' | 
         sort -s -k1,1nr -k2,2nr | 
         cut -f3-
        ) > ~{sub(basename(sv_vcf), "\\.vcf.gz$", "")}.ranked.annotsv.tsv
    >>>

    output {
        File annotsv_annotated_tsv = sub(basename(sv_vcf), "\\.vcf.gz$", "") + ".annotsv.tsv"
        File ranked_annotsv_annotated_tsv = sub(basename(sv_vcf), "\\.vcf.gz$", "") + ".ranked.annotsv.tsv"
    }

    runtime {
        docker: "quay.io/biocontainers/annotsv@sha256:09cc20a86b61fc44b7c1a5d90af8a88fb5e8cacbe56c9938301f2c5fc6ae71fb"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task chord_hrd {
    input {
        File small_variant_vcf
        File sv_vcf
        String pname
        Int threads = 4
    }

    Float file_size = ceil(size(small_variant_vcf, "GB") + size(sv_vcf, "GB"))

    command <<<
    set -euxo pipefail

    # Docker image uses GRIDSS as default. Change to Manta
    sed 's/gridss/manta/g' \
        /opt/chord/extractSigPredictHRD.R > ./extractSigPredictHRD.R
    chmod +x ./extractSigPredictHRD.R

    ./extractSigPredictHRD.R . ~{pname} ~{small_variant_vcf} ~{sv_vcf} 38 2>&1 | tee chord_hrd.log
    rm -f ./extractSigPredictHRD.R
    >>>

    output {
        File chord_log = "chord_hrd.log"
        File chord_prediction = pname + "_chord_prediction.txt"
        File chord_signature = pname + "_chord_signatures.txt"
    }

    runtime {
        docker: "scwatts/chord@sha256:9f6aa44ffefe3f736e66a0e2d7941d4f3e1cc6d848a9a11a17e85a6525e63a77"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}



task prioritize_sv_intogen {
    input {
        File annotSV_tsv
        Int threads
    }

    Float file_size = ceil(size(annotSV_tsv, "GB") + 10)

    command <<<
    set -euxo pipefail
    
    csvtk version

    # Remove any quote from the file
    sed 's/"//g' ~{annotSV_tsv} > ~{basename(annotSV_tsv)}_noquote.tsv

    csvtk join -t \
        ~{basename(annotSV_tsv)}_noquote.tsv \
        /app/Compendium_Cancer_Genes.tsv \
        -f "Gene_name;SYMBOL" |\
            csvtk filter2 -t -f '$Annotation_mode == "split"' |\
            csvtk summary -t -g "$(csvtk headers -t ~{annotSV_tsv} | tr '\n' ',' | sed 's/,$//g')" \
                -f CANCER_TYPE:collapse,COHORT:collapse,TRANSCRIPT:collapse,MUTATIONS:collapse,ROLE:collapse,CGC_GENE:collapse,CGC_CANCER_GENE:collapse,DOMAINS:collapse,2D_CLUSTERS:collapse,3D_CLUSTERS:collapse -s ";" |\
                sed 's/:collapse//g'  > ~{sub(basename(annotSV_tsv), "\\.tsv$", "")}_intogenCCG.tsv

    rm -f ~{basename(annotSV_tsv)}_noquote.tsv

    (head -1 ~{sub(basename(annotSV_tsv), "\\.tsv$", "")}_intogenCCG.tsv | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$5,$6,$14,$26,$54,$81,$104,$106,$107,$118}' && 
     tail -n +2 ~{sub(basename(annotSV_tsv), "\\.tsv$", "")}_intogenCCG.tsv | 
     awk -F '\t' -v OFS='\t' '
     {
        split($107, classifications, /;/);
        max_priority = 0;
        for (i in classifications) {
            gsub(/^[ \t]+|[ \t]+$/, "", classifications[i]);
            if (classifications[i] == "Definitive") curr_priority = 4;
            else if (classifications[i] == "Strong") curr_priority = 3;
            else if (classifications[i] == "Moderate") curr_priority = 2;
            else if (classifications[i] == "Supportive") curr_priority = 1;
            else curr_priority = 0;
            if (curr_priority > max_priority) max_priority = curr_priority;
        }
        print max_priority, ($118+0), $1,$2,$3,$5,$6,$14,$26,$54,$81,$104,$106,$107,$118
     }' | 
     sort -s -k1,1nr -k2,2nr | 
     cut -f3-
    ) > ~{sub(basename(annotSV_tsv), "\\.tsv$", "")}_intogenCCG.ranked.tsv

    >>>

    output {
        File annotSV_intogen_tsv = sub(basename(annotSV_tsv), "\\.tsv$", "") + "_intogenCCG.tsv"
        File annotSV_intogen_ranked_tsv = sub(basename(annotSV_tsv), "\\.tsv$", "") + "_intogenCCG.ranked.tsv"
    }

    runtime {
        docker: "quay.io/pacbio/somatic_general_tools@sha256:a25a2e62b88c73fa3c18a0297654420a4675224eb0cf39fa4192f8a1e92b30d6"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}

task prioritize_small_variants {
    input {
        File vep_annotated_vcf
        Int threads
        String pname="sample"
    }

    Float file_size = ceil(size(vep_annotated_vcf, "GB") + 10)
    String fname = sub(basename(vep_annotated_vcf), "\\.vcf.gz", "") + ".tsv"
    String fname2 = sub(basename(vep_annotated_vcf), "\\.vcf.gz", "") + "_intogenCCG.tsv"

    command <<<
    set -euxo pipefail

    csvtk version

    echo -e "CHROM\tPOS\tREF\tALT\tFORMAT\t~{pname}\t$(bcftools +split-vep ~{vep_annotated_vcf} -l | cut -f2 | tr '\n' '\t' | sed 's/\t$//g')" > ~{fname}
    bcftools +split-vep ~{vep_annotated_vcf} -A tab -f '%CHROM\t%POS\t%REF\t%ALT\t%FORMAT\t%CSQ\n' >> ~{fname}

    csvtk join -t \
        ~{fname} \
        <(sed 's/DOMAINS/CCG_DOMAINS/g' /app/Compendium_Cancer_Genes.tsv) \
        -f SYMBOL |\
            csvtk summary -t -g "$(csvtk headers -t ~{fname} | tr '\n' ',' | sed 's/,$//g')" \
                -f CANCER_TYPE:collapse,COHORT:collapse,TRANSCRIPT:collapse,MUTATIONS:collapse,ROLE:collapse,CGC_GENE:collapse,CGC_CANCER_GENE:collapse,CCG_DOMAINS:collapse,2D_CLUSTERS:collapse,3D_CLUSTERS:collapse -s ";" |\
                sed 's/:collapse//g' > ~{fname2}
    >>>

    output {
        File vep_annotated_tsv = fname
        File vep_annotated_tsv_intogenCCG = fname2
    }

    runtime {
        docker: "quay.io/pacbio/somatic_general_tools@sha256:a25a2e62b88c73fa3c18a0297654420a4675224eb0cf39fa4192f8a1e92b30d6"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
    }
}


