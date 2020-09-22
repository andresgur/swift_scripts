#!/bin/env tclsh
#Fit all the observations of every source in the config file /fit_sources_config_files/sources.txt

# author Andres Gurpide Lasheras, PhD at the institute de recherche en astrophysique et planetologie. email: andres.gurpide@gmail.com
#Created the 25-01-2019

# procedure to read a one column file. Returns the content of the column file as a list
 proc read_file_content {input_file} {
    set input_file [open $input_file r]
    set file_content [read $input_file]
    set list_of_content [split $file_content "\n"]
    close $input_file
    # get rid of last empty line
    set list_of_content [lreplace $list_of_content end end]
    return $list_of_content
    
 }
 
  # create components parameter file output header   
 proc create_header {energy_bands} {
    set header "#params"

    foreach band $energy_bands {
        append header "\t" $band
    }
    append header "\n"
    
    return $header
 }

###Configuration files###
set config_files_dir $::env(HOME)/scripts/swift_scripts/config_model_to_countrate

source ~/scripts/xspec_scripts/xspec_utils.tcl
source ~/scripts/xspec_scripts/loading_utils.tcl
source ~/scripts/xspec_scripts/plot_utils.tcl

chatter 10
query yes

set source_params_file "source_params.txt"

set delimiter "\t"

set swift_file "fit_goodness.log"

set observations_file "observations_info.txt"

set chandra_obs_file "chandra_obs.txt"

set nustar_dir "nustardata"
#chandra data directory
set chandra_dir "chandradata"

set nustar_energy_range "**-3.0 30.0-**"

# energy bands
set energy_bands [list "**-0.3 1.5-**" "**-1.5 10.0-**" "**-0.3 10.0-**"]
set nsimulations 1000

#read arguments
set usage "-m(odel_dir) (model dir from a previous fit run) -s(source) Source to be processed -n Number of simulations to be performed (1000 by default)"

for {set j 0} {$j<$argc} {incr j} {

switch -glob -- [lindex $argv $j] {
-m* {incr j;set model_dir [lindex $argv $j] }
-s* {incr j;set source_dir [lindex $argv $j] }
-n* {incr j;set nsimulations [lindex $argv $j] }
-* {puts "Error unknown option [lindex $argv $j] [lindex $argv [expr $j +1]] $usage"; exit}
}
}
#exit program if parameters were not provided
if {[info exists model_dir]==0 || [info exists source_dir]==0} {
puts "One or more arguments are missing. $usage";exit
}

 tclout chatter
 set chatter_level [lindex $xspec_tclout 0]

# do not echo commands
set xs_echo_script 0

puts "Configuration files: $config_files_dir"

set swift_rmf_arf [read_file_content $config_files_dir/swift_data_paths.txt]
# first is the header
set swift_response [lindex $swift_rmf_arf 1]
set swift_arf [lindex $swift_rmf_arf 2]
set swift_data [lindex $swift_rmf_arf 3]
set swift_back [lindex $swift_rmf_arf 4]
puts "Swift data to be used: \n RMF:$swift_response \n ARF: $swift_arf \n Spectrum: $swift_data"

set TIME_start [clock clicks -milliseconds]
#loop through every source
 #enter source directory
 puts "\n Processing source $source_dir \n"
 if {[file exist $source_dir]} {
  cd $source_dir
  } else {
   puts "ERROR: $source_dir directory not found!"
   exit
  }  
   log "$model_dir/swift_rates.log"
   puts "\n Reading goodness files for source $source_dir... \n"
   set swift_file_contents [read_file $model_dir/$swift_file]

   # iterate each observation
   foreach goodness_info $swift_file_contents {
    
    if {[string match "#*" $goodness_info]} {
     continue
    }
    
    set split_line [split $goodness_info "\t"]
    set epoch [lindex $split_line 0]
    set files [lindex $split_line 1]
    set model_file [lindex $split_line 2]
    set stat chi

    puts "Loading model $model_file and data $files"
    @$model_dir/$files    
    @$model_dir/$model_file
    
    set modelname [regsub -all {_model\.xcm} $model_file {}]
    energies 0.3 50 2000 log
    fit
    #cpd /xs
    # pl eeufs del
    chain unload 1
    parallel walkers 8
    chain walkers 8
    chain length 1000
    chain filetype ascii
    rm chain_$modelname\.ascii
    chain run chain_$modelname\.ascii
    pl eufs
    tclout modpar $modelname
    set npars $xspec_tclout 
    set header "#"
    for {set param 1} {$param <=$npars} {incr param} {
     tclout pinfo $modelname:$param
     scan $xspec_tclout "%s" param_name
     append header "$param_name \t "
    }
    puts $header
    
    set outputtext $header
    puts "Simulating $nsimulations models"
    #Generate parameters sampling from the posterior probability distribution matrix
    for {set i 1} {$i <= $nsimulations} {incr i} {
    tclout simpars
    append outputtext "\n $xspec_tclout"
    puts "$i/$nsimulations"
    }
    #write output
    #file to store posterior distribution from mc sampling
    mkdir $model_dir/swift_models
    set outputdir $model_dir/swift_models
    set outfile $outputdir/simpars_$modelname\_$nsimulations\.dat
    set filemc [open $outfile w]
    puts $filemc $outputtext
    close $filemc
    
    data none
    data $swift_data
    #data none
    resp $env(CALDB)/data/swift/xrt/cpf/rmf/$swift_response
    arf $env(CALDB)/data/swift/xrt/cpf/arf/$swift_arf
    set params_file [open $outfile r]

    tclout modpar $modelname
    set npars $xspec_tclout

    set string_out [create_header $energy_bands]
    
    for {set i 1} {$i <= $nsimulations} {incr i} {
   
    #set parameters out of the mcmh simulations
    set line [gets $params_file]
    if {[string match "#*" $line]} {
         continue
        }
	
   for {set param 1} {$param <=$npars} {incr param} {
    newpar $modelname:$param [lindex $line [expr $param-1]] 
   }
   show model
 # get count rates for each simulation
  append string_out $line 
     foreach band $energy_bands {
       puts "Getting count rate in $band"
       ignore *:$band
       tclout rate 1
       # 2 is the model rate, 0 count rate data, 1 error
       show rate
       set model_rate [format "%.3f" [lindex $xspec_tclout 2] ]
       puts $xspec_tclout
       append string_out "\t" $model_rate
       notice all
       #cpd /xs
       #pl eeufs del
     }
     append string_out "\n"
    }
     set fileout [open "$model_dir/swift_models/swift_rates_$modelname.dat" w]
     puts $fileout $string_out
     close $fileout
     model clear
     data none
    
   }
 cd swift_data/lc
 puts "python ~/scripts/swift_scripts/plot_hr.py --swift_dir ../../$model_dir -t 12"
 python ~/scripts/swift_scripts/plot_hr.py --swift_dir ../../$model_dir -t 12
 



    
 
    
    
 
