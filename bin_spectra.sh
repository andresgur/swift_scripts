# init heasoft
. /home/agurpide/heasoft-6.26.1/x86_64-pc-linux-gnu-libc2.23/headas-init.sh
outbin=XRT_source_opt.fits
for dir in */;do cd $dir 
	back_file=$(echo *pcback.pi)
	source_file=$(echo *pcsource.pi)
	resp_file=$(echo *pc.rmf)
	arf_file=$(echo *pc.arf)
	grppha infile=$source_file outfile=XRT_min20.fits comm='group min 20 & exit' clobber=yes
	fthedit XRT_min20.fits[1] BACKFILE add "'$PWD/$back_file'" comment='Name of background file' longstring = YES chatter=5
	fthedit XRT_min20.fits[1] RESPFILE add "'$PWD/$resp_file'" comment='Name of response file' longstring  = YES chatter=5
	fthedit XRT_min20.fits[1] ANCRFILE add "'$PWD/$arf_file'" comment='Name of ARF file' longstring = YES chatter=5
	ftgrouppha $source_file outfile=$outbin clobber=yes grouptype=opt respfile=$PWD/$resp_file backfile=$PWD/$back_file
	fthedit $outbin BACKFILE add "'$PWD/$back_file'" comment='Name of background file' longstring = YES chatter=5
	fthedit $outbin RESPFILE add "'$PWD/$resp_file'" comment='Name of response file' longstring  = YES chatter=5
	fthedit $outbin ANCRFILE add "'$PWD/$arf_file'" comment='Name of ARF file' longstring = YES chatter=5
	cd .. ;done
