version 1.0 

# Use Severus to call SV
task Severus_sv {
  input {
    File wt_bam
    File wt_bam_index
    File? parental_bam
    File? parental_bam_index
    File trf_bed
    File? phased_vcf
    File? PON_tsv
    String pname = "Sample_1_test"
    Int threads
    Int min_supp_reads
  }

  Float file_size = ceil(size(wt_bam, "GB") + size(parental_bam, "GB") + size(phased_vcf, "GB") + 10)
  command <<<
    set -euxo pipefail
    
    echo "Running Severus for ~{pname}"

    severus --version

    severus \
      --target-bam ~{wt_bam} \
      ~{"--control-bam " + parental_bam} \
      ~{"--phasing-vcf " + phased_vcf} \
      ~{"--PON " + PON_tsv} \
      --out-dir ~{pname + "_severus"} \
      -t ~{threads} \
      --vntr-bed ~{trf_bed} \
      --min-support ~{min_supp_reads} \
      --resolve-overlaps \
      --between-junction-ins \
      --single-bp



    # Compress SVs plots HTML inside somatic_SVs/plots directory
    # Check if the directory exists first
    if [[ -d ~{pname + "_severus/somatic_SVs/plots"} ]]
      then tar -czvf ~{pname + "_severus/somatic_SVs/plots.tar.gz"} ~{pname + "_severus/somatic_SVs/plots"}
    fi
  >>>

  output {
    File? output_vcf = pname + "_severus/somatic_SVs/severus_somatic" + ".vcf"
    File? output_all_vcf = pname + "_severus/all_SVs/severus_all.vcf"
    File? output_breakpoint_clusters = pname + "_severus/somatic_SVs/" + "breakpoint_clusters_list.tsv"
    File? output_breakpoint_clusters_all = pname + "_severus/all_SVs/" + "breakpoint_clusters_list.tsv"
    File? output_somatic_sv_plots = pname + "_severus/somatic_SVs/plots.tar.gz"
  }

  runtime {
    docker: "quay.io/biocontainers/severus@sha256:fb4471e0504d564de78215ae15c081a1bb2022ad51e993eba92bc6fa5052a05d"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

# Use Severus to call SV
task Severus_sv_wt_only {
  input {
    File wt_bam 
    File wt_bam_index
    File trf_bed
    File? phased_vcf
    File? PON_tsv
    String pname = "Sample_1_test"
    Int threads
    Int min_supp_reads
  }

  Float file_size = ceil(size(wt_bam, "GB") + size(phased_vcf, "GB") + 10)

  command <<<
    set -euxo pipefail
    
    echo "Running Severus for ~{pname}"

    severus --version

    severus \
      --target-bam ~{wt_bam} \
      ~{"--phasing-vcf " + phased_vcf} \
      ~{"--PON " + PON_tsv} \
      --out-dir ~{pname + "_severus"} \
      -t ~{threads} \
      --vntr-bed ~{trf_bed} \
      --min-support ~{min_supp_reads} \
      --resolve-overlaps \
      --between-junction-ins \
      --single-bp

    # Compress SVs plots HTML inside somatic_SVs/plots directory
    # Check if the directory exists first
    if [[ -d ~{pname + "_severus/somatic_SVs/plots"} ]]
      then tar -czvf ~{pname + "_severus/somatic_SVs/plots.tar.gz"} ~{pname + "_severus/somatic_SVs/plots"}
    fi
  >>>

  output {
    File output_vcf = pname + "_severus/somatic_SVs/severus_somatic" + ".vcf"
    File output_all_vcf = pname + "_severus/all_SVs/severus_all.vcf"
    File output_breakpoint_clusters = pname + "_severus/somatic_SVs/" + "breakpoint_clusters_list.tsv"
    File output_breakpoint_clusters_all = pname + "_severus/all_SVs/" + "breakpoint_clusters_list.tsv"
    File output_somatic_sv_plots = pname + "_severus/somatic_SVs/plots.tar.gz"
  }

  runtime {
    docker: "quay.io/biocontainers/severus@sha256:fb4471e0504d564de78215ae15c081a1bb2022ad51e993eba92bc6fa5052a05d"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}
task tabix_vcf {
  input {
    File vcf
    File contig_bed
    Int threads
  }

  Float file_size = ceil(size(vcf, "GB") + size(contig_bed, "GB") + 10)
  
  command <<<
    set -euxo pipefail
    echo "indexing ~{vcf}"
    # Get rid of SVLEN=0 and change <INS> to INS
    # to avoid Truvari error problem
    # sed 's/SVLEN=0;//g' ~{vcf} > tmp.vcf
    # sed -i 's/<INS>/INS/g' tmp.vcf

    bcftools --version

    # sort first in case input is not sorted
    bcftools sort \
      -Oz -o tmp.vcf.gz \
      ~{vcf}

    tabix -p vcf tmp.vcf.gz

    # Filter out contigs that are not in the contig bed file
    bcftools view \
      -R ~{contig_bed} \
      -Oz -o ~{basename(vcf) + ".gz"} \
      tmp.vcf.gz

    tabix -p vcf ~{basename(vcf) + ".gz"}
    rm -f tmp.vcf.gz tmp.vcf.gz.tbi
  >>>

  output {
    File output_vcf = basename(vcf) + ".gz"
    File output_vcf_index = basename(vcf) + ".gz.tbi"
  }

  runtime {
    docker: "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

# Use svpack filtering instead of truvari, copied from
# human-wgs-wdl
task svpack_filter_annotated {
	input {
		File sv_vcf

		Array[File] population_vcfs
		Array[File] population_vcf_indices

    Float? svlen

		File gff
	}

	String sv_vcf_basename = basename(sv_vcf, ".vcf.gz")
	Int file_size = ceil(size(sv_vcf, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		echo "svpack version:"
		cat /opt/svpack/.git/HEAD

		svpack \
			filter \
			--pass-only \
			~{"--min-svlen " + svlen} \
			~{sv_vcf} \
			~{sep=' ' prefix('| svpack match -v - ', population_vcfs)} \
		| svpack \
			consequence \
			- \
			~{gff} \
		| svpack \
			tagzygosity \
			- \
		> ~{sv_vcf_basename}.svpack.vcf

		bgzip --version

		bgzip ~{sv_vcf_basename}.svpack.vcf

		tabix --version

		tabix -p vcf ~{sv_vcf_basename}.svpack.vcf.gz
	>>>

	output {
		File output_vcf = "~{sv_vcf_basename}.svpack.vcf.gz"
		File output_vcf_index = "~{sv_vcf_basename}.svpack.vcf.gz.tbi"
	}

	runtime {
		docker: "quay.io/pacbio/svpack@sha256:a680421cb517e1fa4a3097838719a13a6bd655a5e6980ace1b03af9dd707dd75"
		cpu: 4
    memory: "16 GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
	}
}

# Use bcftools to add missing BND mate back after svpack filtering
task recover_mate_bnd {
  input {
    File sv_vcf_original
    File sv_svpack_filtered
  }

  Float file_size = ceil(size(sv_vcf_original, "GB") + size(sv_svpack_filtered, "GB") + 10)

  command <<<
  set -euxo pipefail

  bcftools --version

  # Copy the VCF here and index them
  cp ~{sv_svpack_filtered} ~{basename(sv_svpack_filtered)}
  tabix ~{basename(sv_svpack_filtered)}

  # For original VCF, bgzip if not already, then index
  cp ~{sv_vcf_original} ~{basename(sv_vcf_original)}
  if [[ ~{basename(sv_vcf_original)} == *.vcf ]];
  then
    bgzip ~{sv_vcf_original}
    tabix ~{sv_vcf_original}.gz
  else
    tabix ~{basename(sv_vcf_original)}
  fi

  comm -13 \
    <(bcftools query -f '%ID\t%MATE_ID\n' ~{basename(sv_svpack_filtered)} | rg -v '\.' | cut -f1 | sort) \
    <(bcftools query -f '%ID\t%MATE_ID\n' ~{basename(sv_svpack_filtered)} | rg -v '\.' | cut -f2 | sort) \
    > missing_mate.txt

  # Extract the missing BNDs
  bcftools view \
    -i 'ID=@missing_mate.txt' \
    ~{basename(sv_vcf_original)} |\
      bcftools sort -Oz -o missing_mate.vcf.gz

  # Add the missing BNDs back to the filtered VCF
  bcftools concat \
    ~{basename(sv_svpack_filtered)} \
    missing_mate.vcf.gz |\
      bcftools sort -Oz -o tmp.vcf.gz

  mv tmp.vcf.gz ~{basename(sv_svpack_filtered)}
  rm -f ~{basename(sv_svpack_filtered)}.tbi

  tabix -p vcf ~{basename(sv_svpack_filtered)}
  >>>

  output {
    File output_vcf = basename(sv_svpack_filtered)
    File output_vcf_index = basename(sv_svpack_filtered) + ".tbi"
  }

  runtime {
    docker: "quay.io/biocontainers/bcftools:1.17--h3cc50cf_1"
    cpu: 4
    memory: "16 GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

