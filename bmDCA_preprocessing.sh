#!/bin/bash


while getopts ":wrh" opt; do
  case $opt in
    w)
      # echo "-w: Reiweighting will be applied." 
      REWEIGHTING=1
      ;;
    r)
	  # echo "-r: Converting input from fasta format to numerical format."
	  CONVERT=1
	  ;;
	h)
	  echo " "
	  echo "Preprocessing of data for bmDCA -- Usage: `basename $0` input_alignment -- "
	  echo "Options"
	  echo "  -r: Converts input alignment from fasta format to numerical format, using integers from 0 to 20."
	  echo "  -w: Computes weights of each sequence in alignment. Input should be in numerical format, or -r option used."
	  echo " "
	  exit 0
	  ;;
    \?)
      echo "Invalid option: -$OPTARG" 
      exit 0
      ;;
  esac
done
shift $((OPTIND-1))

mkdir -p Processed
input=$@
out="Processed/msa_numerical.txt"
weights_file="Processed/weights.txt"

##### Converting fasta to numerical alignment
if [ $CONVERT ]; then
	echo " "
	echo "Converting "$input" to numerical format in "$out" ..."
	if test -e $out
	then
		rm $out
		touch $out
	else
		touch $out
	fi
	
	gawk '
	BEGIN{Keep="-ACDEFGHIKLMNPQRSTVWY"; Gap="[BJOUXZ]";printf "">"Processed/temp2";printf "">"Processed/temp1"; Nseq=0; seq_length_tot=0;out="";}
	$0 ~ /^>/ {if(seq_length>0){printf("%s\n",substr(out,1,length(out)-1))>>"Processed/temp1"; if(seq_length>seq_length_tot){seq_length_tot=seq_length;}}
			   out="";seq_length = 0;Nseq++;
			  }
	$0 !~ /^>/ {gsub(Gap,"-",$0);
				for(i=1;i<=length($0);i++){if(index(Keep,substr($0,i,1))){out = out index(Keep,substr($0,i,1))-1 " ";seq_length++;}}
				}
	END {if(seq_length>0){printf("%s\n",out)>>"Processed/temp1";}
		printf("Number of sequences: %d\nSize of sequences: %d\n",Nseq,seq_length_tot);
		printf("%d %d %d\n",Nseq,seq_length,21)>>"Processed/temp2"}
	' $input
	cat Processed/temp2 Processed/temp1 > $out
	rm Processed/temp2 Processed/temp1
	echo "Done!"
	echo " "
else
	cp $1 $out
fi

##### Computing weights
if [ $REWEIGHTING ]
	then
	echo " "
	echo "Computing weights for "$out" ..."
	./sources/reweighting.out $out $weights_file
	echo "Done"
	echo "Output written in "$weights_file"."
	echo " "
fi