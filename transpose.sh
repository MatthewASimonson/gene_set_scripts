# Transpose file row/col:
#########################
# NOTE: this file is space delimited
awk -F ' ' '{
for (f = 1; f <= NF; f++)
a[NR, f] = $f
}
NF > nf { nf = NF }
END {
for (f = 1; f <= nf; f++)
for (r = 1; r <= NR; r++)
printf a[r, f] (r==NR ? RS : FS) 
}' input.transpose  > output.transposed.txt