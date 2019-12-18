#!/bin/bash

CLUSTEREDFASTADIR="/disk2/mmansfield/LCT/Round4/revisions/collection-analysis-stuff/bont-dt-lct/clustered-inputs"

printf "Type\tCluster_Identity\tNum_Seqs\n" > effectors.tb
for USEARCHOUT in $CLUSTEREDFASTADIR/*;
do
	TYPE=$(basename $USEARCHOUT | cut -d '_' -f 1 | cut -d '.' -f 1);
	ID=$(echo $USEARCHOUT | cut -d '.' -f 2,3);
	NUMSEQS=$(grep -c ">" $USEARCHOUT);
	printf "$TYPE\t$ID\t$NUMSEQS\n";
done >> effectors.tb
echo "Done!"
