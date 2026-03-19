#!/bin/bash

getfile() 
{
	mkdir -p /eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/${2}
	scp $T3PSI:/pnfs/psi.ch/cms/trivcat/store/user/mhuwiler/Analysis/submit/${1} /eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/${2}/ #t3psi
}

renameprods()
{
	SUFFIX=${1}
	prods=( D1 D2 D3 D4 D5 A1 A2 A3 A4 A5 B1 B2 B3 B4 B5 )
	for prod in $prods; do
		echo data$prod$SUFFIX.root
	done
}

stripprodname()
{
	SUFFIX=${1}
	for file in $(ls data*${SUFFIX}.root); do
		mv "$file" "`echo $file | sed 's/'${SUFFIX}'//'`"
	done
}

fetchfile() 
{
	scp $T3PSI:/pnfs/psi.ch/cms/trivcat/store/user/mhuwiler/Analysis/submit/hadding/${1} . #t3psi
}
