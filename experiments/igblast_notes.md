* internal files are manually downloaded from here: ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data/human/
  using the following executed from this directory (igblast root)

        $ cat internal_files
        | sed -e 's|^|ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data/human/|'
        | xargs wget -P internal_data/human/
* Get the IMGT databse of V, D, and J genes (imgt/{V,D,J}_imgt.fa) from
  http://www.imgt.org/vquest/refseqh.html
* Clean them up using the edit_imgt_file.pl script taken from
  ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/
  using:

        $ ./edit_imgt_file.pl imgt/imgt_D.fa > imgt/D.fa

  * There are duplicate V genes (for both `IGHV1-45*03` and `IGHV2-70*04`) in
    the output of the above which were resolved by picking the one that agrees
    with:
    https://www.ncbi.nlm.nih.gov/igblast/showGermline.cgi?organism=human&chainType=VH&seqType=nucleotide&functionClass=1

* Install `ncbi-blast+` for `makeblastdb` (or use the one in `bin/` of igblast)
  and create db's for all germline files:

        $ makeblastdb -parse_seqids -dbtype nucl -in imgt/V.fa

  and similar for D and J.
* get mappings via the following from igblast root (`WORKDIR` is my hack)

        $ WORKDIR=../experiments/data/igh-s22 bin/igblastn \
                  -germline_db_V $WORKDIR/imgt/V.fa \
                  -germline_db_D $WORKDIR/imgt/D.fa \
                  -germline_db_J $WORKDIR/imgt/J.fa \
                  -penalty -1 \ # for V ?
                  -D_penalty -3 \
                  -J_penalty -2 \
                  -gapopen -5 \
                  -gapextend -2 \
                  -outfmt 7 \
                  -query $WORKDIR/s22_first_100.fa \
                  -out $WORKDIR/ighblast_first_100.out
  cleanup output via:

        $ OUT=../experiments/data/igh-s22/ighblast_first_100.out \
          grep -e '^[VDJ]' $OUT | \
          awk '{print $1, $2, $3, $4, $5 - $6"/"$5 }' :
          > ../experiments/data/igh-s22/ighblast_first_100_clean.out
