#!/bin/bash

#change these if running in different location
INPUTFASTAS="/disk2/mmansfield/LCT/Round4/revisions/collection-analysis-stuff/bont-dt-lct/inputs"
OUTPUTFASTAS="/disk2/mmansfield/LCT/Round4/revisions/collection-analysis-stuff/bont-dt-lct/clustered-inputs"

for PCTID in 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00;
do
	for PROTEIN in $INPUTFASTAS/*.fa;
	do
		BN=$(basename $PROTEIN | cut -d '.' -f 1);
		echo "    Clustering $BN at $PCTID identity...";
		usearch --cluster_fast $PROTEIN --centroids $OUTPUTFASTAS/$BN.$PCTID --id $PCTID --sizeout > /dev/null 2>&1 ; # clusters every protein at each identity level
	done;
done

echo "Done!"
