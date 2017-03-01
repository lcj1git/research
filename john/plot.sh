#!/bin/bash  
InpFiles='*.gpt'
for gpt in $InpFiles 
do
# is there any input at all 
	if [ -e $gpt ] 
	then
	tmpfile=.tmp-plot.ib3
	errlog=.tmp-err.ib3
# if for some wierd reason you store sth in the funny hidden file then we won't destroy it
		if [ -e $tmpfile -o -e $errlog ] 
		then
			while [[ $yesno != "Y" ]]
			do
				echo "The temp-files exist, overwrite? [Y/n]"
				read yesno
					if [ $yesno = "n" ]
					then
					exit 1
					fi
			done
		fi
# you don't have that file or it is to be overwritten 
# now we read the outputs from the gpt file
	echo `cat $gpt | grep "set output" | awk '{print $3}' | sed s/\"//g` > $tmpfile
# now we can run the plotting with gnuplot 
	gnuplot $gpt
# small test if gnuplot produced the correct output file
	if [ $? -ne 0 ] 
	then
	echo "Gnuplot failed to produce output"
# clean
	rm -f $tmpfile
	exit 1
	fi
# and convert each and single graph into separate files
# first we need to convert eps files to pdfs
		for plot in $(cat $tmpfile) 
		do
		ps2pdf -dEPSCrop ${plot%\.*}"-inc.eps"
# and compile pdfs	
		i=1
                	while [ $i -le 3 ] 
			do
# kill the annoying prompt 
			pdflatex --halt-on-error $plot > $errlog
			if [ $? -ne 0 ] 
# sth went terribly wrong so tail the error log and die 
			then 
# prompt message and clean 			
			echo "pdflatex failed to produce output"
			tail $errlog
		        rm -f ${plot%\.*}{"-inc.eps","-inc.pdf",".log",".aux",".tex"} 
# and clean all the files that have been already produced and not yet compiled 
				for trash in $(cat $tmpfile)
				do
				rm -f ${trash%\.*}{"-inc.eps","-inc.pdf",".tex",".log","aux"}
				done
			rm $errlog $tmpfile
			exit 1
			fi
# or go further 
			i=`echo $i+1 | bc`
			done	
		rm -f ${plot%\.*}{"-inc.eps","-inc.pdf",".log",".aux",".tex"}
		echo "Plot: ${plot%\.*}.pdf created"
		done
# clean this job and take new input 
		rm -f $tmpfile $errlog
	else
	echo "No input files in $PWD directory"
        exit 1
	fi
done
	
exit 0
