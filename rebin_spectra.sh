. $HOME/heasoft-6.26.1/x86_64-pc-linux-gnu-libc2.23/headas-init.sh
spectrum=Medium
spectrum_out=medium_groupped_20.fits
if [[ -d $spectrum ]] ; then cd $spectrum; grppha infile=$spectrum"pc.pi" outfile=$spectrum_out comm='group min 20 & exit' clobber=yes; cd ..;fi
spectrum=High
spectrum_out=high_groupped_20.fits
if [[ -d $spectrum ]]; then cd $spectrum; grppha infile=$spectrum"pc.pi" outfile=$spectrum_out comm='group min 20 & exit' clobber=yes; cd ..;fi
spectrum=Low
spectrum_out=low_groupped_20.fits
if [[ -d $spectrum ]]; then cd $spectrum; grppha infile=$spectrum"pc.pi" outfile=$spectrum_out comm='group min 20 & exit' clobber=yes; cd ..;fi
