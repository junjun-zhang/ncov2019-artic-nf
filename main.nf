#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// include modules
include printHelp from './modules/help.nf'

// import subworkflows
include {articNcovNanopore} from './workflows/articNcovNanoporeExt.nf'
include {ncovIllumina} from './workflows/illuminaNcovExt.nf'
include {ncovIlluminaCram} from './workflows/illuminaNcovExt.nf'
include {songScoreDownload as dnld} from './workflows/song-score-download'
include {songScoreUpload} from './workflows/song-score-upload'
include {getSampleId; prepareFastqPair; performHostFilter; getSecondaryFiles} from './modules/utils'


if (params.help){
    printHelp()
    exit 0
}

if (params.profile){
    println("Profile should have a single dash: -profile")
    System.exit(1)
}

if ( params.illumina ) {
   if ( !params.directory ) {
       println("Please supply a directory containing fastqs or CRAMs with --directory. Specify --cram if supplying a CRAMs directory")
       println("Use --help to print help")
       System.exit(1)
   }
   if ( (params.bed && ! params.ref) || (!params.bed && params.ref) ) {
       println("--bed and --ref must be supplied together")
       System.exit(1)
   }
} else if ( params.nanopolish ) {
   if (! params.basecalled_fastq ) {
       println("Please supply a directory containing basecalled fastqs with --basecalled_fastq. This is the output directory from guppy_barcoder or guppy_basecaller - usually fastq_pass. This can optionally contain barcodeXX directories, which are auto-detected.")
   }
   if (! params.fast5_pass ) {
       println("Please supply a directory containing fast5 files with --fast5_pass (this is the fast5_pass directory)")
   }
   if (! params.sequencing_summary ) {
       println("Please supply the path to the sequencing_summary.txt file from your run with --sequencing_summary")
       System.exit(1)
   }
   if ( params.bed || params.ref ) {
       println("ivarBed and alignerRefPrefix only work in illumina mode")
       System.exit(1)
   }
} else if ( params.medaka ) {
   if (! params.basecalled_fastq ) {
       println("Please supply a directory containing basecalled fastqs with --basecalled_fastq. This is the output directory from guppy_barcoder or guppy_basecaller - usually fastq_pass. This can optionally contain barcodeXX directories, which are auto-detected.")
   }
} else {
       println("Please select a workflow with --nanopolish, --illumina or --medaka, or use --help to print help")
       System.exit(1)
}


if ( ! params.prefix ) {
     println("Please supply a prefix for your output files with --prefix")
     println("Use --help to print help")
     System.exit(1)
} else {
     if ( params.prefix =~ /\// ){
         println("The --prefix that you supplied contains a \"/\", please replace it with another character")
         System.exit(1)
     }
} 



// main workflow
workflow {
   // TODO: download from SONG/SCORE
   dnld(params.study_id, params.analysis_id)
   analysis_metadata = dnld.out.song_analysis
   sequencing_files = dnld.out.files

   getSampleId(analysis_metadata)
   prepareFastqPair(getSampleId.out, sequencing_files)


   if ( params.illumina ) {
       performHostFilter(
          prepareFastqPair.out,
          file(params.human_ref),
          Channel.fromPath(
             getSecondaryFiles(params.human_ref, ['amb', 'ann', 'bwt', 'pac', 'sa']),
             checkIfExists: true
          ).collect()
       )

       if (params.cram) {
           Channel.fromPath( "${params.directory}/**.cram" )
                  .map { file -> tuple(file.baseName, file) }
                  .set{ ch_cramFiles }
       }
       else if (params.analysis_id) {
           ch_filePairs = performHostFilter.out
       }
       else {
	   Channel.fromFilePairs( params.fastqSearchPath, flat: true)
	          .set{ ch_filePairs }
       }
   }
   else {
       // Check to see if we have barcodes
       nanoporeBarcodeDirs = file("${params.basecalled_fastq}/barcode*", type: 'dir', maxdepth: 1 )
       nanoporeNoBarcode = file("${params.basecalled_fastq}/*.fastq", type: 'file', maxdepth: 1)

       if( nanoporeBarcodeDirs ) {
            // Yes, barcodes!
            Channel.fromPath( nanoporeBarcodeDirs )
                   .filter( ~/.*barcode[0-9]{1,4}$/ )
                   .filter{ it.listFiles().size() > 5 }
                   .set{ ch_fastqDirs }
       } else if ( nanoporeNoBarcode ){
            // No, no barcodes
            Channel.fromPath( "${params.basecalled_fastq}", type: 'dir', maxDepth: 1 )
                    .set{ ch_fastqDirs }
      } else {
            println("Couldn't detect whether your Nanopore run was barcoded or not. Use --basecalled_fastq to point to the unmodified guppy output directory.")
            System.exit(1)
      }
   }

   main:
     results = Channel.empty()

     if ( params.nanopolish || params.medaka ) {
         articNcovNanopore(ch_fastqDirs)
         results = articNcovNanopore.out  // TODO: define emits in articNcovNanopore
     } else if ( params.illumina ) {
         if ( params.cram ) {
            ncovIlluminaCram(ch_cramFiles)
            results = ncovIlluminaCram.out  // TODO: define emits in ncovIlluminaCram
         }
         else {
            ncovIllumina(ch_filePairs)
            results = ncovIllumina.out
         }
     } else {
         println("Please select a workflow with --nanopolish, --illumina or --medaka")
     }
     
    // upload results
    // TODO: prepare SONG metadata

    // TODO: upload to SONG/SCORE

   publish:  // this is mainly for testing. All results upload to SONG/SCORE no need for publish
      results.qc_pass to: "${params.outdir}"
      results.sorted_aligned_bam to: "${params.outdir}"
      results.variants to: "${params.outdir}"
      results.consensus_fa to: "${params.outdir}"

}

