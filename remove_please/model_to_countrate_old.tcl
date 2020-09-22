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
    set header "#xmm_obsid\tnustar_obsid\tchandra\tepoch\tmodel\tparams"

    foreach band $energy_bands {
        append header "\t" $band "\t" $band "err" 
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

set fit_goodness_file "fit_goodness.log"

set observations_file "observations_info.txt"

set chandra_obs_file "chandra_obs.txt"

# energy bands
set energy_bands [list "**-0.3 1.5-**" "**-1.5 10.0-**" "**-0.3 10.0-**"]

#read arguments
set usage "-m(odel_dir) (model dir from a previous fit run) -d(ata) Some fiducial data so XSPEC does not crash when trying to load a RMF and ARF alone -s(source) Source to be processed"

for {set j 0} {$j<$argc} {incr j} {

switch -glob -- [lindex $argv $j] {
-m* {incr j;set model_dir [lindex $argv $j] }
-s* {incr j;set source_dir [lindex $argv $j] }
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

    #get nustar and xmm observations files
    puts "\n Reading goodness files for source $source_dir... \n"
    set goodness_contents [read_file $model_dir/$fit_goodness_file]
    set string_out [create_header $energy_bands]
    # iterate each observation
    foreach goodness_info $goodness_contents {
     
     if {[string match "#*" $goodness_info]} {
      continue
     }
     
     set split_line [split $goodness_info "\t"]
     set epoch [lindex $split_line 0]
     set model [lindex $split_line 8]
     puts "Model found $model"
     
     data none
     model none
     # set model
     @$model
     # load any data to fool XSPEC
     data $swift_data
     #data none
     resp $env(CALDB)/data/swift/xrt/cpf/rmf/$swift_response
     arf $env(CALDB)/data/swift/xrt/cpf/arf/$swift_arf
     
     append string_out $xmm_obs_id "\t" $nustar_obs_id "\t" $chandra_obs_id "\t" $epoch "\t" $model "\t" $params 
     foreach band $energy_bands {
       puts "Getting count rate in $band"
       ignore *:$band
       tclout rate 1
       # 2 is the model rate, 0 count rate data, 1 error
       show rate
       set model_rate [format "%.3f" [lindex $xspec_tclout 2] ]
       puts $xspec_tclout
       append string_out "\t" $model_rate "\t" [expr {sqrt($model_rate)}]
       cpd /xs
       pl eufs del
       notice all
       
       #cpd /xs
       #pl eeufs del
     }
     model none 
     data none
     append string_out "\n"
     
 }
     set fileout [open "$model_dir/swift_rates.dat" w]
     puts $fileout $string_out
     close $fileout

log none
puts "Results stored into folder $model_dir"
set TIME_taken [expr [clock clicks -milliseconds] - $TIME_start]
puts "Time taken [expr $TIME_taken/1000/60/60] h"


    
 
    
    
 
