#!/bin/bash

chromopainter=../bin/chromopainter
data=../example

${chromopainter} \
	-g ${data}/BrahuiYorubaSimulationChrom22.haplotypes \
	-r ${data}/BrahuiYorubaSimulationChrom22.recomrates \
	-t ${data}/BrahuiYorubaSimulation.idfile.txt \
	-f ${data}/BrahuiYorubaSimulation.poplist.txt 0 0 \
	-o ${data}/BrahuiYorubaSimulationChrom22 \
	-b 
