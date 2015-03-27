##########################################################################################
# Use UNIX commands to extract individual ID from cel .list files and arrange into columns
##########################################################################################

# Clear any previously created files:
rm IIDs.rough; rm IIDs.final;rm Plates.rough; rm Plate.names.final; rm Plate.rough; rm Plate.list.final

ls *.list | while read x; do cut -d "_" -f 11 $x >> IIDs.rough ; done # insert IID with extension
cut -d "." -f 1 IIDs.rough >> IIDs.final # remove extensions and leave just IIDs

ls *.list | while read x; do cut -d "/" -f 8 $x >> Plates.rough ; done # insert plate name with extension
cut -d "_" -f 1 Plates.rough >> Plate.names.final # remove extensions and leave just IIDs # remove extensions and leave just plate names

# Now paste together two files into single 2 column file:

paste IIDs.final Plate.names.final > Plate.rough
grep -v cel Plate.rough > Plate.list.final
