version 1.1

import "humanwgs_structs.wdl"
import "wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration
import "upstream/upstream.wdl" as Upstream
import "joint/joint.wdl" as Joint
import "downstream/downstream.wdl" as Downstream
import "wdl-common/wdl/tasks/bcftools.wdl" as Bcftools
import "wdl-common/wdl/tasks/trgt.wdl" as Trgt
import "wdl-common/wdl/tasks/write_ped_phrank.wdl" as Write_ped_phrank
import "tertiary/tertiary.wdl" as TertiaryAnalysis
import "wdl-common/wdl/tasks/utilities.wdl" as Utilities
import "somatic_ports/somatic_annotation.wdl" as Somatic_annotation
import "somatic_ports/somatic_calling.wdl" as Somatic_calling
workflow humanwgs_family {
  meta {
    description: "PacBio HiFi human whole genome sequencing pipeline, with joint calling for related samples."
  }

  parameter_meta {
    family: {
      name: "Family struct describing samples, relationships, and unaligned BAM paths"
    }
    ref_map_file: {
      name: "TSV containing reference genome file paths; must match backend"
    }
    deepvariant_version: {
      name: "DeepVariant version"
    }
    custom_deepvariant_model_tar: {
      name: "Custom DeepVariant model tarball"
    }
    pharmcat_version: {
      name: "PharmCAT version"
    }
    pharmcat_min_coverage: {
      name: "Minimum coverage for PharmCAT"
    }
    phenotypes: {
      name: "Comma-delimited list of HPO codes for phenotypes"
    }
    tertiary_map_file: {
      name: "TSV containing tertiary analysis file paths and thresholds; must match backend"
    }
    glnexus_mem_gb: {
      name: "Override GLnexus memory request (GB)"
    }
    pbsv_call_mem_gb: {
      name: "Override PBSV call memory request (GB)"
    }
    gpu: {
      name: "Use GPU when possible"
    }
    backend: {
      name: "Backend where the workflow will be executed",
      choices: ["GCP", "Azure", "AWS-HealthOmics", "HPC"]
    }
    zones: {
      name: "Zones where compute will take place; required if backend is set to 'GCP'"
    }
    gpuType: {
      name: "GPU type to use; required if gpu is set to `true` for cloud backends; must match backend"
    }
    container_registry: {
      name: "Container registry where workflow images are hosted. If left blank, PacBio's public Quay.io registry will be used. Must be set if backend is set to 'AWS-HealthOmics'"
    }
    preemptible: {
      name: "Where possible, run tasks preemptibly"
    }
    debug_version: {
      name: "Debug version for testing purposes"
    }
  }

  input {
    Family family

    File ref_map_file

    # These options are only intended for testing purposes.
    # There is no guarantee that the pipeline will work with
    # other version of DeepVariant or with custom models.
    String deepvariant_version = "1.6.1"
    File? custom_deepvariant_model_tar

    String pharmcat_version = "2.15.4"
    Int pharmcat_min_coverage = 10

    String phenotypes = "HP:0000001"
    File? tertiary_map_file

    Int? glnexus_mem_gb
    Int? pbsv_call_mem_gb

    Boolean gpu = false

    # Backend configuration
    String backend
    String? zones
    String? gpuType
    String? container_registry

    Boolean preemptible = true

    String? debug_version
    # Somatic Ports
    ## SV-calling
    File trf_bed = "/hpc/scratch/group.data.science/ram_temp/HiFi-human-WGS-editing-QC-WDL/human_GRCh38_no_alt_analysis_set.trf.bed"
    File ref_bed = "/hpc/scratch/group.data.science/ram_temp/HiFi-human-WGS-editing-QC-WDL/chr.bed"
    File ref_gff = "/hpc/scratch/group.data.science/ram_temp/HiFi-human-WGS-editing-QC-WDL/Homo_sapiens.GRCh38.112.chr.reformatted.gff3"
    File control_vcf = "/hpc/scratch/group.data.science/ram_temp/HiFi-human-WGS-editing-QC-WDL/severus.jasmine.AN10.AC4.nosample.vcf.gz"
    File control_vcf_index = "/hpc/scratch/group.data.science/ram_temp/HiFi-human-WGS-editing-QC-WDL/severus.jasmine.AN10.AC4.nosample.vcf.gz.tbi"
    File severus_pon_tsv = "/hpc/scratch/group.data.science/ram_temp/HiFi-human-WGS-editing-QC-WDL/PoN_1000G_hg38_extended.tsv.gz"
    ## Annotation 
    File vep_cache = "/hpc/scratch/group.data.science/ram_temp/HiFi-human-WGS-editing-QC-WDL/homo_sapiens_refseq_vep_112_GRCh38.tar.gz"
    File annotsv_cache = "/hpc/scratch/group.data.science/ram_temp/HiFi-human-WGS-editing-QC-WDL/annotsv_cache.tar.gz"
  }

  call BackendConfiguration.backend_configuration {
    input:
      backend            = backend,
      zones              = zones,
      gpuType            = gpuType,
      container_registry = container_registry
  }

  RuntimeAttributes default_runtime_attributes = if preemptible then backend_configuration.spot_runtime_attributes else backend_configuration.on_demand_runtime_attributes

  Map [String, String] ref_map = read_map(ref_map_file)

  Boolean single_sample = length(family.samples) == 1

  #############################################################
  # 1) Per‐sample Upstream & Downstream for Wild‐Type (WT)     #
  #############################################################





  ####################################
  # 1a) UPSTREAM: per‐sample calls  #
  ####################################
  scatter (sample in family.samples) {
      String sid_wt = sample.sample_id      
    ## Wild-type branch (always run)
      call Upstream.upstream as upstream_wt {
            input:
              sample_id                    = "~{sid_wt}_WT",
              sex                          = sample.sex,
              hifi_reads                   = sample.wild_type_bams,
              ref_map_file                 = ref_map_file,
              deepvariant_version          = deepvariant_version,
              custom_deepvariant_model_tar = custom_deepvariant_model_tar,
              single_sample                = single_sample,
              gpu                          = gpu,
              default_runtime_attributes   = default_runtime_attributes
      }
      #call Bcftools.bcftools_merge_assembly_align as merge_sv_vcfs_align_assembly {
      #          input:
      #              vcfs = [upstream_wt.large_sv_filtered_vcf, select_first([upstream_wt.sv_vcf])],
      #              out_prefix  = "~{sid_wt}.merged_structural_variants",
      #              runtime_attributes = default_runtime_attributes
      #}

  } 
    ####################################
    # 1b)         JOINT:                #
    ####################################
    # after your upstream scatter and single_sample Boolean…
  if (!single_sample) {
          call Joint.joint as joint_wild {
              input:
                  family_id                  = family.family_id,
                  sample_ids                 = sid_wt,
                  gvcfs                      = upstream_wt.small_variant_gvcf,
                  gvcf_indices               = upstream_wt.small_variant_gvcf_index,
                  svsigs                     = flatten(upstream_wt.svsigs),
                  ref_map_file               = ref_map_file,
                  glnexus_mem_gb             = glnexus_mem_gb,
                  pbsv_call_mem_gb           = pbsv_call_mem_gb,
                  default_runtime_attributes = default_runtime_attributes
          }
  }
  
  scatter (sample_index in range(length(family.samples))) { 
      call Bcftools.bcftools_merge_assembly_align as merge_sv_vcfs_align_assembly {
                input:
                    vcfs = [upstream_wt.large_sv_filtered_vcf[sample_index], select_all(select_first([joint_wild.split_joint_structural_variant_vcfs,upstream_wt.sv_vcf]))[sample_index]],
                    out_prefix  = "~{sid_wt[sample_index]}_WT.merged_structural_variants",
                    runtime_attributes = default_runtime_attributes
      }
  }  
  
    ####################################
    # 1c) DOWNSTREAM: per‐sample calls  #
    ####################################
  scatter (sample_index in range(length(family.samples))) { 
      call Downstream.downstream as downstream_wt {
          input:
          sample_id                  = "~{sid_wt[sample_index]}_WT",
          small_variant_vcf          = select_first([joint_wild.split_joint_small_variant_vcfs, upstream_wt.small_variant_vcf])[sample_index],
          small_variant_vcf_index    = select_first([joint_wild.split_joint_small_variant_vcf_indices, upstream_wt.small_variant_vcf_index])[sample_index],
          sv_vcf                     = merge_sv_vcfs_align_assembly.merged_vcf[sample_index],
          sv_vcf_index               = merge_sv_vcfs_align_assembly.merged_vcf_index[sample_index],
          trgt_vcf                   = upstream_wt.trgt_vcf[sample_index],
          trgt_vcf_index             = upstream_wt.trgt_vcf_index[sample_index],
          aligned_bam                = upstream_wt.out_bam[sample_index],
          aligned_bam_index          = upstream_wt.out_bam_index[sample_index],
          pharmcat_version           = pharmcat_version,
          pharmcat_min_coverage      = pharmcat_min_coverage,
          ref_map_file               = ref_map_file,
          default_runtime_attributes = default_runtime_attributes
      }
  }
  Map[String, Array[String]] stats_wt = {
          'sample_id': sid_wt,
          'num_reads': upstream_wt.stat_num_reads,
          'read_length_mean': upstream_wt.stat_read_length_mean,
          'read_length_median': upstream_wt.stat_read_length_median,
          'read_quality_mean': upstream_wt.stat_read_quality_mean,
          'read_quality_median': upstream_wt.stat_read_quality_median,
          'mapped_read_count': downstream_wt.stat_mapped_read_count,
          'mapped_percent': downstream_wt.stat_mapped_percent,
          'mean_depth': upstream_wt.stat_mean_depth,
          'inferred_sex': upstream_wt.inferred_sex,
          'stat_phased_basepairs': downstream_wt.stat_phased_basepairs,
          'phase_block_ng50': downstream_wt.stat_phase_block_ng50,
          'cpg_combined_count': downstream_wt.stat_combined_cpg_count,
          'cpg_hap1_count': downstream_wt.stat_hap1_cpg_count,
          'cpg_hap2_count': downstream_wt.stat_hap2_cpg_count,
          'SNV_count': downstream_wt.stat_SNV_count,
          'TSTV_ratio': downstream_wt.stat_TSTV_ratio,
          'HETHOM_ratio': downstream_wt.stat_HETHOM_ratio,
          'INDEL_count': downstream_wt.stat_INDEL_count,
          'sv_DUP_count': downstream_wt.stat_sv_DUP_count,
          'sv_DEL_count': downstream_wt.stat_sv_DEL_count,
          'sv_INS_count': downstream_wt.stat_sv_INS_count,
          'sv_INV_count': downstream_wt.stat_sv_INV_count,
          'sv_BND_count': downstream_wt.stat_sv_BND_count,
          'cnv_DUP_count': upstream_wt.stat_cnv_DUP_count,
          'cnv_DEL_count': upstream_wt.stat_cnv_DEL_count,
          'cnv_DUP_sum': upstream_wt.stat_cnv_DUP_sum,
          'cnv_DEL_sum': upstream_wt.stat_cnv_DEL_sum,
          'trgt_genotyped_count': upstream_wt.stat_trgt_genotyped_count,
          'trgt_uncalled_count': upstream_wt.stat_trgt_uncalled_count
  }
  
  call Utilities.consolidate_stats as consolidate_stats_wt {
       input:
           id                 = family.family_id,
           stats              = stats_wt,
           runtime_attributes = default_runtime_attributes
  }



  #############################################################
  # 2) Per‐sample Parental (only run if parental_bams exist)  #
  #############################################################
  # A one‐per‐sample scatter that produces an Array[File?] of IDs
  #   none if no parental_bams, or "<sample_id>_PARENT" otherwise
  scatter (s in family.samples) {
    Sample? maybe_parent_samples = 
        if (length(s.parental_bams) > 0) then
            s
        else
            None
  }   

  
  # Now filter out the None values to get only the parental sample IDs
  Array[Sample] parental_samples = select_all(maybe_parent_samples)
  if (length(parental_samples) > 0) {
    ####################################
    # 2a) UPSTREAM: per‐sample calls  #
    ####################################
    scatter (sample in parental_samples) {
        String sid_parental = sample.sample_id      
        if (length(sample.parental_bams) > 0) {
            call Upstream.upstream as upstream_parental {
                input:
                  sample_id                    = "~{sid_parental}_PARENT",
                  sex                          = sample.sex,
                  hifi_reads                   = sample.parental_bams,
                  ref_map_file                 = ref_map_file,
                  deepvariant_version          = deepvariant_version,
                  custom_deepvariant_model_tar = custom_deepvariant_model_tar,
                  single_sample                = single_sample,
                  gpu                          = gpu,
                  default_runtime_attributes   = default_runtime_attributes
            }
            
        }
            
    }
    

    ####################################
    # 2b)         JOINT:                #
    ####################################
    # after your upstream scatter and single_sample Boolean…
          
    if (!single_sample) {
      call Joint.joint as joint_parental {
        input:
          family_id                  = family.family_id,
          sample_ids                 = sid_parental,
          gvcfs                      = select_all(upstream_parental.small_variant_gvcf),
          gvcf_indices               = select_all(upstream_parental.small_variant_gvcf_index),
          svsigs                     = flatten(select_all(upstream_parental.svsigs)),
          ref_map_file               = ref_map_file,
          glnexus_mem_gb             = glnexus_mem_gb,
          pbsv_call_mem_gb           = pbsv_call_mem_gb,
          default_runtime_attributes = default_runtime_attributes
      }
    }
    scatter (sample_index in range(length(parental_samples))) { 
          call Bcftools.bcftools_merge_assembly_align as merge_sv_vcfs_align_assembly_parental {
                    input:
                        vcfs = [select_all(upstream_parental.large_sv_filtered_vcf)[sample_index], select_all(select_first([joint_parental.split_joint_structural_variant_vcfs,upstream_parental.sv_vcf]))[sample_index]],
                        out_prefix  = "~{sid_parental[sample_index]}_PARENT.merged_structural_variants",
                        runtime_attributes = default_runtime_attributes
          }
      }  

    ####################################
    # 2c) DOWNSTREAM: per‐sample calls  #
    ####################################
    # Use struct for better organization
    # Resolve arrays based on single_sample flag 
    scatter (sample_index in range(length(parental_samples))) {
        if (length(parental_samples[sample_index].parental_bams) > 0) {
            Array[File] resolved_joint_vcfs = select_first([joint_parental.split_joint_small_variant_vcfs, []])
            Array[File] resolved_upstream_vcfs = select_all(upstream_parental.small_variant_vcf) 
            Array[File] resolved_joint_vcfs_index = select_first([joint_parental.split_joint_small_variant_vcf_indices, []])
            Array[File] resolved_upstream_vcfs_index = select_all(upstream_parental.small_variant_vcf_index)

        
            Array[File] final_vcfs = if (!single_sample) then resolved_joint_vcfs else resolved_upstream_vcfs
            Array[File] final_vcf_indices = if (!single_sample) then resolved_joint_vcfs_index else resolved_upstream_vcfs_index
            
            #Array[File] final_vcfs = select_first([resolved_joint_vcfs, resolved_upstream_vcfs])
            #Array[File] final_vcf_indices = select_first([resolved_joint_vcfs_index, resolved_upstream_vcfs_index])
            
            call Downstream.downstream as downstream_parental {
                input:
                    sample_id                  = "~{sid_parental[sample_index]}_PARENT",
                    small_variant_vcf = final_vcfs[sample_index], 
                    small_variant_vcf_index = final_vcf_indices[sample_index], 
                    sv_vcf                     = select_all(merge_sv_vcfs_align_assembly_parental.merged_vcf)[sample_index],
                    sv_vcf_index               = select_all(merge_sv_vcfs_align_assembly_parental.merged_vcf_index)[sample_index],
                    trgt_vcf                   = select_first([upstream_parental.trgt_vcf[sample_index]]),
                    trgt_vcf_index             = select_first([upstream_parental.trgt_vcf_index[sample_index]]),
                    aligned_bam                = select_first([upstream_parental.out_bam[sample_index]]),
                    aligned_bam_index          = select_first([upstream_parental.out_bam_index[sample_index]]),
                    pharmcat_version           = pharmcat_version,
                    pharmcat_min_coverage      = pharmcat_min_coverage,
                    ref_map_file               = ref_map_file,
                    default_runtime_attributes = default_runtime_attributes
            }
        }
        
    }
 
    Map[String, Array[String]] stats_parental = {
      'sample_id': select_all(select_first([sid_parental, []])),
      'num_reads': select_all(select_first([upstream_parental.stat_num_reads, []])),
      'read_length_mean': select_all(select_first([upstream_parental.stat_read_length_mean, []])),
      'read_length_median': select_all(select_first([upstream_parental.stat_read_length_median, []])),
      'read_quality_mean': select_all(select_first([upstream_parental.stat_read_quality_mean, []])),
      'read_quality_median': select_all(select_first([upstream_parental.stat_read_quality_median, []])),
      'mapped_read_count': select_all(select_first([downstream_parental.stat_mapped_read_count, []])),
      'mapped_percent': select_all(select_first([downstream_parental.stat_mapped_percent, []])),
      'mean_depth': select_all(select_first([upstream_parental.stat_mean_depth, []])),
      'inferred_sex': select_all(select_first([upstream_parental.inferred_sex, []])),
      'stat_phased_basepairs': select_all(select_first([downstream_parental.stat_phased_basepairs, []])),
      'phase_block_ng50': select_all(select_first([downstream_parental.stat_phase_block_ng50, []])),
      'cpg_combined_count': select_all(select_first([downstream_parental.stat_combined_cpg_count, []])),
      'cpg_hap1_count': select_all(select_first([downstream_parental.stat_hap1_cpg_count, []])),
      'cpg_hap2_count': select_all(select_first([downstream_parental.stat_hap2_cpg_count, []])),
      'SNV_count': select_all(select_first([downstream_parental.stat_SNV_count, []])),
      'TSTV_ratio': select_all(select_first([downstream_parental.stat_TSTV_ratio, []])),
      'HETHOM_ratio': select_all(select_first([downstream_parental.stat_HETHOM_ratio, []])),
      'INDEL_count': select_all(select_first([downstream_parental.stat_INDEL_count, []])),
      'sv_DUP_count': select_all(select_first([downstream_parental.stat_sv_DUP_count, []])),
      'sv_DEL_count': select_all(select_first([downstream_parental.stat_sv_DEL_count, []])),
      'sv_INS_count': select_all(select_first([downstream_parental.stat_sv_INS_count, []])),
      'sv_INV_count': select_all(select_first([downstream_parental.stat_sv_INV_count, []])),
      'sv_BND_count': select_all(select_first([downstream_parental.stat_sv_BND_count, []])),
      'cnv_DUP_count': select_all(select_first([upstream_parental.stat_cnv_DUP_count, []])),
      'cnv_DEL_count': select_all(select_first([upstream_parental.stat_cnv_DEL_count, []])),
      'cnv_DUP_sum': select_all(select_first([upstream_parental.stat_cnv_DUP_sum, []])),
      'cnv_DEL_sum': select_all(select_first([upstream_parental.stat_cnv_DEL_sum, []])),
      'trgt_genotyped_count': select_all(select_first([upstream_parental.stat_trgt_genotyped_count, []])),
      'trgt_uncalled_count': select_all(select_first([upstream_parental.stat_trgt_uncalled_count, []]))
    }
    call Utilities.consolidate_stats as consolidate_stats_parental {
      input:
          id                 = family.family_id,
          stats              = stats_parental,
          runtime_attributes = default_runtime_attributes
    }
  }

    #############################################################
    # 3)                  Tertiary Analysis                     #
    #############################################################

    #############################################################
    #       MERGING SVs and small variants PER SAMPLE           #
    #############################################################



  if (length(parental_samples) > 0) {
        scatter (sample_index_wt in range(length(family.samples))) {
            scatter (sample_index_parental in range(length(parental_samples))) {
                if (sid_wt[sample_index_wt] == sid_parental[sample_index_parental]) {
                    call Bcftools.bcftools_merge as merge_small_variant_vcfs_per_sample {
                           input:
                             vcfs = select_all(flatten([[downstream_wt.phased_small_variant_vcf[sample_index_wt]], if defined(downstream_parental.phased_small_variant_vcf) then [downstream_parental.phased_small_variant_vcf[sample_index_parental]] else []])),
                             vcf_indices        = select_all(flatten([[downstream_wt.phased_small_variant_vcf_index[sample_index_wt]], if defined(downstream_parental.phased_small_variant_vcf_index) then [downstream_parental.phased_small_variant_vcf_index[sample_index_parental]] else []])),
                             out_prefix         = "~{sid_wt[sample_index_wt]}.joint.~{ref_map['name']}.small_variants.phased",
                             runtime_attributes = default_runtime_attributes
                    }

                    call Bcftools.bcftools_merge as merge_sv_vcfs_per_sample {
                       input:
                             vcfs = select_all(flatten([[downstream_wt.phased_sv_vcf[sample_index_wt]], if defined(downstream_parental.phased_sv_vcf) then [downstream_parental.phased_sv_vcf[sample_index_parental]] else []])),
                         vcf_indices        = select_all(flatten([[downstream_wt.phased_sv_vcf_index[sample_index_wt]], if defined(downstream_parental.phased_sv_vcf) then [downstream_parental.phased_sv_vcf_index[sample_index_parental]] else []])),
                         out_prefix         = "~{sid_wt[sample_index_wt]}.joint.~{ref_map['name']}.structural_variants.phased",
                         runtime_attributes = default_runtime_attributes
                    }
                }
            }
        }
  }   
  

    #############################################################
    #       MERGING SVs and small variants ACROSS SAMPLES       #
    #############################################################



  if (!single_sample) {
     call Bcftools.bcftools_merge as merge_small_variant_vcfs {
       input:
         vcfs               = flatten([downstream_wt.phased_small_variant_vcf, select_all(select_first([downstream_parental.phased_small_variant_vcf, []]))]),
         vcf_indices        = flatten([downstream_wt.phased_small_variant_vcf_index,select_all(select_first([downstream_parental.phased_small_variant_vcf_index, []]))]),
         out_prefix         = "~{family.family_id}.joint.~{ref_map['name']}.small_variants.phased",
         runtime_attributes = default_runtime_attributes
     }

     call Bcftools.bcftools_merge as merge_sv_vcfs {
       input:
         vcfs               = flatten([downstream_wt.phased_sv_vcf, select_all(select_first([downstream_parental.phased_sv_vcf, []]))]),
         vcf_indices        = flatten([downstream_wt.phased_sv_vcf_index,select_all(select_first([downstream_parental.phased_sv_vcf_index, []]))]),
         out_prefix         = "~{family.family_id}.joint.~{ref_map['name']}.structural_variants.phased",
         runtime_attributes = default_runtime_attributes
     }

     call Trgt.trgt_merge {
       input:
         vcfs               = flatten([downstream_wt.phased_trgt_vcf, select_all(select_first([downstream_parental.phased_trgt_vcf, []]))]),
         vcf_indices        = flatten([downstream_wt.phased_trgt_vcf_index,select_all(select_first([downstream_parental.phased_trgt_vcf_index, []]))]),
         ref_fasta          = ref_map["fasta"],                              # !FileCoercion
         ref_index          = ref_map["fasta_index"],                        # !FileCoercion
         out_prefix         = "~{family.family_id}.merged.~{ref_map['name']}.trgt",
         runtime_attributes = default_runtime_attributes
     }
  }  




  #############################################################
  #                  SEVERUS & ANNOTATION                     #
  #############################################################


  if (defined(tertiary_map_file)) {
    
    scatter (sample in family.samples) {
      Array[File] hifi_reads = flatten([sample.parental_bams,sample.wild_type_bams])
    }

    call Write_ped_phrank.write_ped_phrank {
      input:
        id                 = family.family_id,
        family             = family,
        phenotypes         = phenotypes,
        disk_size          = ceil(size(flatten(hifi_reads), "GB")) + 10,
        runtime_attributes = default_runtime_attributes
    }

    ####################################
    # 3a)         SEVERUS RUN (WT):    #
    ####################################
    scatter (sample_index in range(length(family.samples))) {
        call Somatic_calling.Severus_sv_wt_only as phased_severus_wt {
            input:
                wt_bam        = downstream_wt.merged_haplotagged_bam[sample_index],
                wt_bam_index  = downstream_wt.merged_haplotagged_bam_index[sample_index],
                trf_bed          = trf_bed,
                phased_vcf       = downstream_wt.phased_small_variant_vcf[sample_index],
                threads          = 32,
                min_supp_reads   = 3,
                PON_tsv          = severus_pon_tsv
        }


        call Somatic_calling.tabix_vcf as tabix_vcf_wt {
         input:
           vcf        = phased_severus_wt.output_vcf,
           contig_bed = ref_bed,
           threads    = 32
        } 
        call Somatic_calling.svpack_filter_annotated as svpack_filter_annotated_wt {
         input:
           sv_vcf               = tabix_vcf_wt.output_vcf,
           population_vcfs      = [control_vcf],
           population_vcf_indices = [control_vcf_index],
           gff                  = ref_gff
        }
        call Somatic_calling.recover_mate_bnd as recover_mate_bnd_wt {
         input:
           sv_vcf_original    = phased_severus_wt.output_vcf,
           sv_svpack_filtered = svpack_filter_annotated_wt.output_vcf
        }


        call Somatic_annotation.vep_annotate as annotateGermline {
         input:
            input_vcf           = select_first([merge_small_variant_vcfs.merged_vcf, downstream_wt.phased_small_variant_vcf[sample_index]]),
            vep_cache           = select_first([vep_cache]),
            ref_fasta           = ref_map["fasta"],
            ref_fasta_index     = ref_map["fasta_index"],
            threads             = 32
        }


        call Somatic_annotation.annotsv as annotateSV {
            input: 
                sv_vcf              = select_first([merge_sv_vcfs.merged_vcf, downstream_wt.phased_sv_vcf[sample_index]]),
                sv_vcf_index        = select_first([merge_sv_vcfs.merged_vcf_index, downstream_wt.phased_sv_vcf_index[sample_index]]),
                annotsv_cache       = select_first([annotsv_cache]),  
                threads             = 32
        }


        call Somatic_annotation.annotsv as annotateSeverusSVfiltered {
            input: 
                sv_vcf              = select_first([recover_mate_bnd_wt.output_vcf]),
                sv_vcf_index        = select_first([recover_mate_bnd_wt.output_vcf_index ]),
                annotsv_cache       = select_first([annotsv_cache]),  
                threads             = 32
        }

        call Somatic_annotation.prioritize_sv_intogen as prioritize_sv {
            input: 
                annotSV_tsv             = annotateSV.annotsv_annotated_tsv,
                threads                 = 32
        } 

        
        call Somatic_annotation.prioritize_sv_intogen as prioritize_Severus {
            input: 
                annotSV_tsv             = annotateSeverusSVfiltered.annotsv_annotated_tsv, 
                threads                 = 32
        } 

        #call Somatic_annotation.prioritize_small_variants as prioritizeSomatic {
        #    input: 
        #        vep_annotated_vcf       = annotateGermline.vep_annotated_vcf[sample_index_wt], 
        #        threads                 = 32
        #}





    }


    
    ####################################
    # 3b)   SEVERUS RUN (PARENTAL):    #
    ####################################
    
    if (length(parental_samples) > 0) {
        scatter (sample_index_wt in range(length(family.samples))) {
            scatter (sample_index_parental in range(length(parental_samples))) {
                if (sid_wt[sample_index_wt] == sid_parental[sample_index_parental]) {
                    call Somatic_calling.Severus_sv as phased_severus_pair {
                        input:
                            wt_bam = downstream_wt.merged_haplotagged_bam[sample_index_wt],
                            wt_bam_index  = downstream_wt.merged_haplotagged_bam_index[sample_index_wt],
                            parental_bam     = downstream_parental.merged_haplotagged_bam[sample_index_parental],
                            parental_bam_index = downstream_parental.merged_haplotagged_bam_index[sample_index_parental],
                            trf_bed          = trf_bed,
                            phased_vcf       = downstream_wt.phased_small_variant_vcf[sample_index_wt],
                            threads          = 32,
                            min_supp_reads   = 3,
                            PON_tsv          = severus_pon_tsv
                    }


                    call Somatic_calling.tabix_vcf as tabix_vcf_pair {
                     input:
                       vcf        = select_first([phased_severus_pair.output_vcf]),
                       contig_bed = ref_bed,
                       threads    = 32
                    } 
                    call Somatic_calling.svpack_filter_annotated as svpack_filter_annotated_pair {
                    input:
                       sv_vcf               = tabix_vcf_pair.output_vcf,
                       population_vcfs      = [control_vcf],
                       population_vcf_indices = [control_vcf_index],
                       gff                  = ref_gff
                    }
                    call Somatic_calling.recover_mate_bnd as recover_mate_bnd_pair {
                    input:
                       sv_vcf_original    = select_first([phased_severus_pair.output_vcf]),
                       sv_svpack_filtered = svpack_filter_annotated_pair.output_vcf
                    }
                    

                    Array[File] small_variant_vcfs = if (defined(merge_small_variant_vcfs_per_sample.merged_vcf)) then select_all(flatten(select_first([merge_small_variant_vcfs_per_sample.merged_vcf, []]))) else downstream_wt.phased_small_variant_vcf
                    
                    call Somatic_annotation.vep_annotate as annotateGermline_pair {
                     input:
                        input_vcf           = select_first([merge_small_variant_vcfs.merged_vcf, small_variant_vcfs[sample_index_wt]]),
                        vep_cache           = select_first([vep_cache]),
                        ref_fasta           = ref_map["fasta"],
                        ref_fasta_index     = ref_map["fasta_index"],
                        threads             = 32
                    }
                    

                    Array[File] sv_vcfs = if (defined(merge_sv_vcfs_per_sample.merged_vcf)) then select_all(flatten(select_first([merge_sv_vcfs_per_sample.merged_vcf, []]))) else downstream_wt.phased_sv_vcf
                    Array[File] sv_vcfs_index = if (defined(merge_sv_vcfs_per_sample.merged_vcf_index)) then select_all(flatten(select_first([merge_sv_vcfs_per_sample.merged_vcf_index, []]))) else downstream_wt.phased_sv_vcf_index

                    call Somatic_annotation.annotsv as annotateSV_pair {
                        input: 
                            sv_vcf              = select_first([merge_sv_vcfs.merged_vcf, sv_vcfs[sample_index_wt]]),
                            sv_vcf_index        = select_first([merge_sv_vcfs.merged_vcf_index, sv_vcfs_index[sample_index_wt]]),
                            annotsv_cache       = select_first([annotsv_cache]),  
                            threads             = 32
                    }


                    call Somatic_annotation.annotsv as annotateSeverusSVfiltered_pair {
                        input: 
                            sv_vcf              = select_first([recover_mate_bnd_pair.output_vcf]),
                            sv_vcf_index        = select_first([recover_mate_bnd_pair.output_vcf_index]),
                            annotsv_cache       = select_first([annotsv_cache]),  
                            threads             = 32
                    }

                    call Somatic_annotation.prioritize_sv_intogen as prioritize_sv_pair {
                        input: 
                            annotSV_tsv             = annotateSV_pair.annotsv_annotated_tsv,
                            threads                 = 32
                    } 

                    
                    call Somatic_annotation.prioritize_sv_intogen as prioritize_Severus_pair {
                        input: 
                            annotSV_tsv             = annotateSeverusSVfiltered_pair.annotsv_annotated_tsv, 
                            threads                 = 32
                    } 

                   # call Somatic_annotation.prioritize_small_variants as prioritizeSomatic_pair {
                    #    input: 
                    #        vep_annotated_vcf       = annotateGermline.vep_annotated_vcf[sample_index_wt], 
                    #        threads                 = 32
                    #}

                } 
            }
        }
        
    }
    
    ####################################
    # 3c)    Tertiary Analysis RUN     #
    ####################################
 

#   call TertiaryAnalysis.tertiary_analysis {
#        input:
#            pedigree                   = write_ped_phrank.pedigree,
#            phrank_lookup              = write_ped_phrank.phrank_lookup,
#            small_variant_vcf          = select_first([merge_small_variant_vcfs.merged_vcf, downstream_wt.phased_small_variant_vcf[0]]),
#            small_variant_vcf_index    = select_first([merge_small_variant_vcfs.merged_vcf_index, downstream_wt.phased_small_variant_vcf_index[1]]),
#            sv_vcf                     = select_first([merge_sv_vcfs.merged_vcf, downstream_wt.phased_sv_vcf[0]]),
#            sv_vcf_index               = select_first([merge_sv_vcfs.merged_vcf_index, downstream_wt.phased_sv_vcf_index[0]]),
#            ref_map_file               = ref_map_file,
#            tertiary_map_file          = select_first([tertiary_map_file]),
#            default_runtime_attributes = default_runtime_attributes
#    }
#  }
#



  }


  output {
# to maintain order of samples
    Array[String] sample_ids = flatten([sid_wt,select_all(select_first([sid_parental, []]))])
    File stats_file_wt          = consolidate_stats_wt.output_tsv
    File? stats_file_parental          = consolidate_stats_parental.output_tsv
    
    ## WT OUTPUTS
    # bam stats
    Array[File]   bam_stats_wt                = upstream_wt.read_length_and_quality
    Array[File]   read_length_plot_wt         = upstream_wt.read_length_plot
    Array[File?]  read_quality_plot_wt        = upstream_wt.read_quality_plot
    Array[String] stat_num_reads_wt           = upstream_wt.stat_num_reads
    Array[String] stat_read_length_mean_wt    = upstream_wt.stat_read_length_mean
    Array[String] stat_read_length_median_wt  = upstream_wt.stat_read_length_median
    Array[String] stat_read_quality_mean_wt   = upstream_wt.stat_read_quality_mean
    Array[String] stat_read_quality_median_wt = upstream_wt.stat_read_quality_median

    # merged, haplotagged alignments
    Array[File]   merged_haplotagged_bam_wt       = downstream_wt.merged_haplotagged_bam
    Array[File]   merged_haplotagged_bam_index_wt = downstream_wt.merged_haplotagged_bam_index
    Array[String] stat_mapped_read_count_wt       = downstream_wt.stat_mapped_read_count
    Array[String] stat_mapped_percent_wt          = downstream_wt.stat_mapped_percent
    Array[File]   mapq_distribution_plot_wt       = downstream_wt.mapq_distribution_plot
    Array[File]   mg_distribution_plot_wt         = downstream_wt.mg_distribution_plot

    # mosdepth outputs
    Array[File]   mosdepth_summary_wt                 = upstream_wt.mosdepth_summary
    Array[File]   mosdepth_region_bed_wt              = upstream_wt.mosdepth_region_bed
    Array[File]   mosdepth_region_bed_index_wt        = upstream_wt.mosdepth_region_bed_index
    Array[File]   mosdepth_depth_distribution_plot_wt = upstream_wt.mosdepth_depth_distribution_plot
    Array[String] stat_mean_depth_wt                  = upstream_wt.stat_mean_depth
    Array[String] inferred_sex_wt                     = upstream_wt.inferred_sex

    #Array[File]  merged_assembly_aligned_sv_vcf       = upstream_wt.merged_assembly_aligned_sv_vcf

    #Array[File]  merged_assembly_aligned_sv_vcf_index       = upstream_wt.merged_assembly_aligned_sv_vcf_index
    # phasing stats
    Array[File]   phase_stats_wt           = downstream_wt.phase_stats
    Array[File]   phase_blocks_wt          = downstream_wt.phase_blocks
    Array[File]   phase_haplotags_wt       = downstream_wt.phase_haplotags
    Array[String] stat_phased_basepairs_wt = downstream_wt.stat_phased_basepairs
    Array[String] stat_phase_block_ng50_wt = downstream_wt.stat_phase_block_ng50

    # cpg_pileup outputs
    Array[File?]  cpg_combined_bed_wt        = downstream_wt.cpg_combined_bed
    Array[File?]  cpg_combined_bed_index_wt  = downstream_wt.cpg_combined_bed_index
    Array[File?]  cpg_hap1_bed_wt            = downstream_wt.cpg_hap1_bed
    Array[File?]  cpg_hap1_bed_index_wt      = downstream_wt.cpg_hap1_bed_index
    Array[File?]  cpg_hap2_bed_wt            = downstream_wt.cpg_hap2_bed
    Array[File?]  cpg_hap2_bed_index_wt      = downstream_wt.cpg_hap2_bed_index
    Array[File?]  cpg_combined_bw_wt         = downstream_wt.cpg_combined_bw
    Array[File?]  cpg_hap1_bw_wt             = downstream_wt.cpg_hap1_bw
    Array[File?]  cpg_hap2_bw_wt             = downstream_wt.cpg_hap2_bw
    Array[String] stat_cpg_hap1_count_wt     = downstream_wt.stat_hap1_cpg_count
    Array[String] stat_cpg_hap2_count_wt     = downstream_wt.stat_hap2_cpg_count
    Array[String] stat_cpg_combined_count_wt = downstream_wt.stat_combined_cpg_count

    # sv outputs
    Array[File] phased_sv_vcf_wt       = downstream_wt.phased_sv_vcf
    Array[File] phased_sv_vcf_index_wt = downstream_wt.phased_sv_vcf_index

    # sv stats
    Array[String] stat_sv_DUP_count_wt = downstream_wt.stat_sv_DUP_count
    Array[String] stat_sv_DEL_count_wt = downstream_wt.stat_sv_DEL_count
    Array[String] stat_sv_INS_count_wt = downstream_wt.stat_sv_INS_count
    Array[String] stat_sv_INV_count_wt = downstream_wt.stat_sv_INV_count
    Array[String] stat_sv_BND_count_wt = downstream_wt.stat_sv_BND_count

    # small variant outputs
    Array[File] phased_small_variant_vcf_wt       = downstream_wt.phased_small_variant_vcf
    Array[File] phased_small_variant_vcf_index_wt = downstream_wt.phased_small_variant_vcf_index
    Array[File] small_variant_gvcf_wt             = upstream_wt.small_variant_gvcf
    Array[File] small_variant_gvcf_index_wt       = upstream_wt.small_variant_gvcf_index

    # small variant stats
    Array[File]   small_variant_stats_wt             = downstream_wt.small_variant_stats
    Array[File]   bcftools_roh_out_wt                = downstream_wt.bcftools_roh_out
    Array[File]   bcftools_roh_bed_wt                = downstream_wt.bcftools_roh_bed
    Array[String] stat_small_variant_SNV_count_wt    = downstream_wt.stat_SNV_count
    Array[String] stat_small_variant_INDEL_count_wt  = downstream_wt.stat_INDEL_count
    Array[String] stat_small_variant_TSTV_ratio_wt   = downstream_wt.stat_TSTV_ratio
    Array[String] stat_small_variant_HETHOM_ratio_wt = downstream_wt.stat_HETHOM_ratio
    Array[File]   snv_distribution_plot_wt           = downstream_wt.snv_distribution_plot
    Array[File]   indel_distribution_plot_wt         = downstream_wt.indel_distribution_plot

    # trgt outputs
    Array[File]   phased_trgt_vcf_wt           = downstream_wt.phased_trgt_vcf
    Array[File]   phased_trgt_vcf_index_wt     = downstream_wt.phased_trgt_vcf_index
    Array[File]   trgt_spanning_reads_wt       = upstream_wt.trgt_spanning_reads
    Array[File]   trgt_spanning_reads_index_wt = upstream_wt.trgt_spanning_reads_index
    Array[File]   trgt_coverage_dropouts_wt    = upstream_wt.trgt_coverage_dropouts
    Array[String] stat_trgt_genotyped_count_wt = upstream_wt.stat_trgt_genotyped_count
    Array[String] stat_trgt_uncalled_count_wt  = upstream_wt.stat_trgt_uncalled_count

    # paraphase outputs
    Array[File]  paraphase_output_json_wt         = upstream_wt.paraphase_output_json
    Array[File]  paraphase_realigned_bam_wt       = upstream_wt.paraphase_realigned_bam
    Array[File]  paraphase_realigned_bam_index_wt = upstream_wt.paraphase_realigned_bam_index
    Array[File?] paraphase_vcfs_wt                = upstream_wt.paraphase_vcfs

    # per sample cnv outputs
    Array[File]   cnv_vcf_wt              = upstream_wt.cnv_vcf
    Array[File]   cnv_vcf_index_wt        = upstream_wt.cnv_vcf_index
    Array[File]   cnv_copynum_bedgraph_wt = upstream_wt.cnv_copynum_bedgraph
    Array[File]   cnv_depth_bw_wt         = upstream_wt.cnv_depth_bw
    Array[File]   cnv_maf_bw_wt           = upstream_wt.cnv_maf_bw
    Array[String] stat_cnv_DUP_count_wt   = upstream_wt.stat_cnv_DUP_count
    Array[String] stat_cnv_DEL_count_wt   = upstream_wt.stat_cnv_DEL_count
    Array[String] stat_cnv_DUP_sum_wt     = upstream_wt.stat_cnv_DUP_sum
    Array[String] stat_cnv_DEL_sum_wt     = upstream_wt.stat_cnv_DEL_sum

    # PGx outputs
    Array[File]  pbstarphase_json_wt        = downstream_wt.pbstarphase_json
    Array[File?] pharmcat_match_json_wt     = downstream_wt.pharmcat_match_json
    Array[File?] pharmcat_phenotype_json_wt = downstream_wt.pharmcat_phenotype_json
    Array[File?] pharmcat_report_html_wt    = downstream_wt.pharmcat_report_html
    Array[File?] pharmcat_report_json_wt    = downstream_wt.pharmcat_report_json

    # Assembly outputs
    #bam to fasta outputs
    Array[File] fasta_output_wt = upstream_wt.fasta_output
    
    
    # hifiasm upstream_wt outputs
    Array[File] asm_1_wt = upstream_wt.asm_1
    Array[File] asm_2_wt = upstream_wt.asm_2

    # pav outputs
    Array[File] pav_vcf_wt = upstream_wt.pav_vcf
    Array[File] pav_vcf_index_wt = upstream_wt.pav_vcf_index

    #liftover outputs
    Array[File] lifted_vcf_wt = upstream_wt.lifted_vcf
    Array[File] reject_vcf_wt = upstream_wt.reject_vcf 

    #LARGE SV FILTER OUTPUTS
    Array[File] large_sv_filtered_vcf_wt = upstream_wt.large_sv_filtered_vcf

    Array[File] large_sv_filtered_vcf_index_wt = upstream_wt.large_sv_filtered_vcf_index
    Array[File] merged_assembly_aligned_sv_vcf_wt = select_first([merge_sv_vcfs_align_assembly.merged_vcf])
    Array[File] merged_assembly_aligned_sv_vcf_index_wt = select_first([merge_sv_vcfs_align_assembly.merged_vcf_index])    
    
    ## PARENTAL OUTPUTS
    # bam stats
    Array[File]   bam_stats_parental                = select_all(select_first([upstream_parental.read_length_and_quality, []]))
    Array[File]   read_length_plot_parental         = select_all(select_first([upstream_parental.read_length_plot, []]))
    Array[File?]  read_quality_plot_parental        = select_all(select_first([upstream_parental.read_quality_plot, []]))
    Array[String] stat_num_reads_parental           = select_all(select_first([upstream_parental.stat_num_reads, []]))
    Array[String] stat_read_length_mean_parental    = select_all(select_first([upstream_parental.stat_read_length_mean, []]))
    Array[String] stat_read_length_median_parental  = select_all(select_first([upstream_parental.stat_read_length_median, []]))
    Array[String] stat_read_quality_mean_parental   = select_all(select_first([upstream_parental.stat_read_quality_mean, []]))
    Array[String] stat_read_quality_median_parental = select_all(select_first([upstream_parental.stat_read_quality_median, []]))

    # merged, haplotagged alignments
    Array[File]   merged_haplotagged_bam_parental       = select_all(select_first([downstream_parental.merged_haplotagged_bam, []]))
    Array[File]   merged_haplotagged_bam_index_parental = select_all(select_first([downstream_parental.merged_haplotagged_bam_index, []]))
    Array[String] stat_mapped_read_count_parental       = select_all(select_first([downstream_parental.stat_mapped_read_count, []]))
    Array[String] stat_mapped_percent_parental          = select_all(select_first([downstream_parental.stat_mapped_percent, []]))
    Array[File]   mapq_distribution_plot_parental       = select_all(select_first([downstream_parental.mapq_distribution_plot, []]))
    Array[File]   mg_distribution_plot_parental         = select_all(select_first([downstream_parental.mg_distribution_plot, []]))

    # mosdepth outputs
    Array[File]   mosdepth_summary_parental                 = select_all(select_first([upstream_parental.mosdepth_summary, []]))
    Array[File]   mosdepth_region_bed_parental              = select_all(select_first([upstream_parental.mosdepth_region_bed, []]))
    Array[File]   mosdepth_region_bed_index_parental        = select_all(select_first([upstream_parental.mosdepth_region_bed_index, []]))
    Array[File]   mosdepth_depth_distribution_plot_parental = select_all(select_first([upstream_parental.mosdepth_depth_distribution_plot, []]))
    Array[String] stat_mean_depth_parental                  = select_all(select_first([upstream_parental.stat_mean_depth, []]))
    Array[String] inferred_sex_parental                     = select_all(select_first([upstream_parental.inferred_sex, []]))

    # phasing stats
    Array[File]   phase_stats_parental           = select_all(select_first([downstream_parental.phase_stats, []]))
    Array[File]   phase_blocks_parental          = select_all(select_first([downstream_parental.phase_blocks, []]))
    Array[File]   phase_haplotags_parental       = select_all(select_first([downstream_parental.phase_haplotags, []]))
    Array[String] stat_phased_basepairs_parental = select_all(select_first([downstream_parental.stat_phased_basepairs, []]))
    Array[String] stat_phase_block_ng50_parental = select_all(select_first([downstream_parental.stat_phase_block_ng50, []]))

    # cpg_pileup outputs
    Array[File?]  cpg_combined_bed_parental        = select_all(select_first([downstream_parental.cpg_combined_bed, []]))
    Array[File?]  cpg_combined_bed_index_parental  = select_all(select_first([downstream_parental.cpg_combined_bed_index, []]))
    Array[File?]  cpg_hap1_bed_parental            = select_all(select_first([downstream_parental.cpg_hap1_bed, []]))
    Array[File?]  cpg_hap1_bed_index_parental      = select_all(select_first([downstream_parental.cpg_hap1_bed_index, []]))
    Array[File?]  cpg_hap2_bed_parental            = select_all(select_first([downstream_parental.cpg_hap2_bed, []]))
    Array[File?]  cpg_hap2_bed_index_parental      = select_all(select_first([downstream_parental.cpg_hap2_bed_index, []]))
    Array[File?]  cpg_combined_bw_parental         = select_all(select_first([downstream_parental.cpg_combined_bw, []]))
    Array[File?]  cpg_hap1_bw_parental             = select_all(select_first([downstream_parental.cpg_hap1_bw, []]))
    Array[File?]  cpg_hap2_bw_parental             = select_all(select_first([downstream_parental.cpg_hap2_bw, []]))
    Array[String] stat_cpg_hap1_count_parental     = select_all(select_first([downstream_parental.stat_hap1_cpg_count, []]))
    Array[String] stat_cpg_hap2_count_parental     = select_all(select_first([downstream_parental.stat_hap2_cpg_count, []]))
    Array[String] stat_cpg_combined_count_parental = select_all(select_first([downstream_parental.stat_combined_cpg_count, []]))

    # sv outputs
    Array[File] phased_sv_vcf_parental       = select_all(select_first([downstream_parental.phased_sv_vcf, []]))
    Array[File] phased_sv_vcf_index_parental = select_all(select_first([downstream_parental.phased_sv_vcf_index, []]))

    # sv stats
    Array[String] stat_sv_DUP_count_parental = select_all(select_first([downstream_parental.stat_sv_DUP_count, []]))
    Array[String] stat_sv_DEL_count_parental = select_all(select_first([downstream_parental.stat_sv_DEL_count, []]))
    Array[String] stat_sv_INS_count_parental = select_all(select_first([downstream_parental.stat_sv_INS_count, []]))
    Array[String] stat_sv_INV_count_parental = select_all(select_first([downstream_parental.stat_sv_INV_count, []]))
    Array[String] stat_sv_BND_count_parental = select_all(select_first([downstream_parental.stat_sv_BND_count, []]))

    # small variant outputs
    Array[File] phased_small_variant_vcf_parental       = select_all(select_first([downstream_parental.phased_small_variant_vcf, []]))
    Array[File] phased_small_variant_vcf_index_parental = select_all(select_first([downstream_parental.phased_small_variant_vcf_index, []]))
    Array[File] small_variant_gvcf_parental             = select_all(select_first([upstream_parental.small_variant_gvcf, []]))
    Array[File] small_variant_gvcf_index_parental       = select_all(select_first([upstream_parental.small_variant_gvcf_index, []]))

    # small variant stats
    Array[File]   small_variant_stats_parental             = select_all(select_first([downstream_parental.small_variant_stats, []]))
    Array[File]   bcftools_roh_out_parental                = select_all(select_first([downstream_parental.bcftools_roh_out, []]))
    Array[File]   bcftools_roh_bed_parental                = select_all(select_first([downstream_parental.bcftools_roh_bed, []]))
    Array[String] stat_small_variant_SNV_count_parental    = select_all(select_first([downstream_parental.stat_SNV_count, []]))
    Array[String] stat_small_variant_INDEL_count_parental  = select_all(select_first([downstream_parental.stat_INDEL_count, []]))
    Array[String] stat_small_variant_TSTV_ratio_parental   = select_all(select_first([downstream_parental.stat_TSTV_ratio, []]))
    Array[String] stat_small_variant_HETHOM_ratio_parental = select_all(select_first([downstream_parental.stat_HETHOM_ratio, []]))
    Array[File]   snv_distribution_plot_parental           = select_all(select_first([downstream_parental.snv_distribution_plot, []]))
    Array[File]   indel_distribution_plot_parental         = select_all(select_first([downstream_parental.indel_distribution_plot, []]))

    # trgt outputs
    Array[File]   phased_trgt_vcf_parental           = select_all(select_first([downstream_parental.phased_trgt_vcf, []]))
    Array[File]   phased_trgt_vcf_index_parental     = select_all(select_first([downstream_parental.phased_trgt_vcf_index, []]))
    Array[File]   trgt_spanning_reads_parental       = select_all(select_first([upstream_parental.trgt_spanning_reads, []]))
    Array[File]   trgt_spanning_reads_index_parental = select_all(select_first([upstream_parental.trgt_spanning_reads_index, []]))
    Array[File]   trgt_coverage_dropouts_parental    = select_all(select_first([upstream_parental.trgt_coverage_dropouts, []]))
    Array[String] stat_trgt_genotyped_count_parental = select_all(select_first([upstream_parental.stat_trgt_genotyped_count, []]))
    Array[String] stat_trgt_uncalled_count_parental  = select_all(select_first([upstream_parental.stat_trgt_uncalled_count, []]))

    # paraphase outputs
    Array[File]  paraphase_output_json_parental         = select_all(select_first([upstream_parental.paraphase_output_json, []]))
    Array[File]  paraphase_realigned_bam_parental       = select_all(select_first([upstream_parental.paraphase_realigned_bam, []]))
    Array[File]  paraphase_realigned_bam_index_parental = select_all(select_first([upstream_parental.paraphase_realigned_bam_index, []]))
    Array[File?] paraphase_vcfs_parental                = select_all(select_first([upstream_parental.paraphase_vcfs, []]))

    # per sample cnv outputs
    Array[File]   cnv_vcf_parental              = select_all(select_first([upstream_parental.cnv_vcf, []]))
    Array[File]   cnv_vcf_index_parental        = select_all(select_first([upstream_parental.cnv_vcf_index, []]))
    Array[File]   cnv_copynum_bedgraph_parental = select_all(select_first([upstream_parental.cnv_copynum_bedgraph, []]))
    Array[File]   cnv_depth_bw_parental         = select_all(select_first([upstream_parental.cnv_depth_bw, []]))
    Array[File]   cnv_maf_bw_parental           = select_all(select_first([upstream_parental.cnv_maf_bw, []]))
    Array[String] stat_cnv_DUP_count_parental   = select_all(select_first([upstream_parental.stat_cnv_DUP_count, []]))
    Array[String] stat_cnv_DEL_count_parental   = select_all(select_first([upstream_parental.stat_cnv_DEL_count, []]))
    Array[String] stat_cnv_DUP_sum_parental     = select_all(select_first([upstream_parental.stat_cnv_DUP_sum, []]))
    Array[String] stat_cnv_DEL_sum_parental     = select_all(select_first([upstream_parental.stat_cnv_DEL_sum, []]))

    # PGx outputs
    Array[File]  pbstarphase_json_parental        = select_all(select_first([downstream_parental.pbstarphase_json, []]))
    Array[File?] pharmcat_match_json_parental     = select_all(select_first([downstream_parental.pharmcat_match_json, []]))
    Array[File?] pharmcat_phenotype_json_parental = select_all(select_first([downstream_parental.pharmcat_phenotype_json, []]))
    Array[File?] pharmcat_report_html_parental    = select_all(select_first([downstream_parental.pharmcat_report_html, []]))
    Array[File?] pharmcat_report_json_parental    = select_all(select_first([downstream_parental.pharmcat_report_json, []]))

    # MERGED WT+PARENTAL call outputs PER SAMPLE
    Array[File?] parental_wt_small_variants_vcf_per_sample       = select_all(flatten(select_first([merge_small_variant_vcfs_per_sample.merged_vcf, []])))
    Array[File?] parental_wt_small_variants_vcf_index = select_all(flatten(select_first([merge_small_variant_vcfs_per_sample.merged_vcf_index, []])))
    Array[File?] parental_wt_sv_vcf                   = select_all(flatten(select_first([merge_sv_vcfs_per_sample.merged_vcf, []])))
    Array[File?] parental_wt_sv_vcf_index             = select_all(flatten(select_first([merge_sv_vcfs_per_sample.merged_vcf_index, []])))
    
    # hifiasm upstream_parental outputs
    Array[File?] asm_1_parental = select_all(select_first([upstream_parental.asm_1]))

    Array[File?] asm_2_parental = select_all(select_first([upstream_parental.asm_2]))

    # pav outputs
    Array[File?] pav_vcf_parental = select_all(select_first([upstream_parental.pav_vcf]))
    Array[File?] pav_vcf_index_parental = select_all(select_first([upstream_parental.pav_vcf_index]))

    #liftover outputs
    Array[File?] lifted_vcf_parental = select_all(select_first([upstream_parental.lifted_vcf]))
    Array[File?] reject_vcf_parental = select_all(select_first([upstream_parental.reject_vcf ]))

    #LARGE SV FILTER OUTPUTS
    Array[File?] large_sv_filtered_vcf_parental = select_all(select_first([upstream_parental.large_sv_filtered_vcf]))
    Array[File?] large_sv_filtered_vcf_index_parental = select_all(select_first([upstream_parental.large_sv_filtered_vcf_index]))
    
    Array[File?] merged_assembly_aligned_sv_vcf_parental = select_all(select_first([merge_sv_vcfs_align_assembly_parental.merged_vcf]))
    Array[File?] merged_assembly_aligned_sv_vcf_index_parental = select_all(select_first([merge_sv_vcfs_align_assembly_parental.merged_vcf_index]))    

    # joint call outputs
    File? joint_small_variants_vcf       = merge_small_variant_vcfs.merged_vcf
    File? joint_small_variants_vcf_index = merge_small_variant_vcfs.merged_vcf_index
    File? joint_sv_vcf                   = merge_sv_vcfs.merged_vcf
    File? joint_sv_vcf_index             = merge_sv_vcfs.merged_vcf_index
    File? joint_trgt_vcf                 = trgt_merge.merged_vcf
    File? joint_trgt_vcf_index           = trgt_merge.merged_vcf_index

 # tertiary analysis outputs
    ## WT OUTPUTS
    File? pedigree                                      = write_ped_phrank.pedigree
 #    File? tertiary_small_variant_filtered_vcf           = tertiary_analysis.small_variant_filtered_vcf
 #    File? tertiary_small_variant_filtered_vcf_index     = tertiary_analysis.small_variant_filtered_vcf_index
 #    File? tertiary_small_variant_filtered_tsv           = tertiary_analysis.small_variant_filtered_tsv
 #    File? tertiary_small_variant_compound_het_vcf       = tertiary_analysis.small_variant_compound_het_vcf
 #    File? tertiary_small_variant_compound_het_vcf_index = tertiary_analysis.small_variant_compound_het_vcf_index
 #    File? tertiary_small_variant_compound_het_tsv       = tertiary_analysis.small_variant_compound_het_tsv
 #    File? tertiary_sv_filtered_vcf                      = tertiary_analysis.sv_filtered_vcf
 #    File? tertiary_sv_filtered_vcf_index                = tertiary_analysis.sv_filtered_vcf_index
 #    File? tertiary_sv_filtered_tsv                      = tertiary_analysis.sv_filtered_tsv
    
#    #Somatic SV calling 
    Array[File] Severus_somatic_vcf_wt                           = select_first([phased_severus_wt.output_vcf])
    Array[File] Severus_all_vcf_wt                               = select_first([phased_severus_wt.output_all_vcf])
    Array[File] Severus_breakpoint_cluster_wt                    = select_first([phased_severus_wt.output_breakpoint_clusters])
    Array[File] Severus_breakpoint_cluster_all_wt                = select_first([phased_severus_wt.output_breakpoint_clusters_all])
    Array[File] Severus_cluster_plots_wt                         = select_first([phased_severus_wt.output_somatic_sv_plots])
#
    Array[File] Severus_tabix_vcf_wt                             = select_first([tabix_vcf_wt.output_vcf])
    Array[File] Severus_filtered_vcf_wt = select_first([recover_mate_bnd_wt.output_vcf])
    Array[File] Severus_filtered_vcf_index_wt = select_first([recover_mate_bnd_wt.output_vcf_index])
  


    #Somatic annotation analysis outputs
    Array[File] vep_annotated_vcf                             = select_first([annotateGermline.vep_annotated_vcf])
    Array[File] annotsv_annotated_tsv                         = select_first([annotateSV.annotsv_annotated_tsv])
    Array[File] Severus_annotated_tsv                         = select_first([annotateSeverusSVfiltered.annotsv_annotated_tsv])
    Array[File] annotsv_annotated_tsv_intogen                         = select_first([prioritize_sv.annotSV_intogen_tsv])
    Array[File] Severus_annotated_tsv_intogen                         = select_first([prioritize_Severus.annotSV_intogen_tsv])
#    
#
#    File? small_variant_tsv_annotated                           = prioritizeSomatic.vep_annotated_tsv
#    File? small_variant_tsv_CCG                                 = prioritizeSomatic.vep_annotated_tsv_intogenCCG


    ## PAREENTAL OUTPUTS
    Array[File?] Severus_somatic_vcf_pair                           = select_all(flatten(select_first([phased_severus_pair.output_vcf, []])))
    Array[File?] Severus_all_vcf_pair                               = select_all(flatten(select_first([phased_severus_pair.output_all_vcf, []])))
    Array[File?] Severus_breakpoint_cluster_pair                    = select_all(flatten(select_first([phased_severus_pair.output_breakpoint_clusters, []])))
    Array[File?] Severus_breakpoint_cluster_all_pair                = select_all(flatten(select_first([phased_severus_pair.output_breakpoint_clusters_all, []])))
    Array[File?] Severus_cluster_plots_pair                         = select_all(flatten(select_first([phased_severus_pair.output_somatic_sv_plots, []])))
#
    Array[File?] Severus_tabix_vcf_pair                             = select_all(flatten(select_first([tabix_vcf_pair.output_vcf, []])))
    Array[File?] Severus_filtered_vcf_pair = select_all(flatten(select_first([recover_mate_bnd_pair.output_vcf, []])))
    Array[File?] Severus_filtered_vcf_index_pair = select_all(flatten(select_first([recover_mate_bnd_pair.output_vcf_index, []])))
  


    #Somatic annotation analysis outputs
    Array[File?] vep_annotated_vcf_pair                             = select_all(flatten(select_first([annotateGermline_pair.vep_annotated_vcf, []])))
    Array[File?] annotsv_annotated_tsv_pair                         = select_all(flatten(select_first([annotateSV_pair.annotsv_annotated_tsv, []])))
    Array[File?] Severus_annotated_tsv_pair                         = select_all(flatten(select_first([annotateSeverusSVfiltered_pair.annotsv_annotated_tsv, []])))
    Array[File?] annotsv_annotated_tsv_intogen_pair                         = select_all(flatten(select_first([prioritize_sv_pair.annotSV_intogen_tsv, []])))
    Array[File?] Severus_annotated_tsv_intogen_pair                         = select_all(flatten(select_first([prioritize_Severus_pair.annotSV_intogen_tsv, []])))
#    
#
#    File_pair? small_variant_tsv_annotated                           _pair= prioritizeSomatic.vep_annotated_tsv
#    File? small_variant_tsv_CCG                                 = prioritizeSomatic.vep_annotated_tsv_intogenCCG

    # workflow metadata
    String workflow_name    = "humanwgs_family"
    String workflow_version = "v2.0.7" + if defined(debug_version) then "~{"-" + debug_version}" else ""
  }
}
