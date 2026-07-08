#!/usr/bin/env bash

read -p "Enter Observation ID: " OBSID

echo "OBSID: $OBSID"

echo "Running Resolve Pipeline..."
# Run the pipeline
xapipeline indir="$OBSID" outdir="$OBSID"_rsl_reproc steminputs=xa"$OBSID" stemoutputs=DEFAULT entry_stage=1 exit_stage=3 instrument=resolve verify_input=no create_ehkmkf=no calc_pointing=yes calc_optaxis=yes calc_gtilost=no calc_adrgti=no calc_mxsgti=no rsl_gainfile=CALDB calmethod=FE55 linetocorrect=MnKa seed=650504 numevent=1000 minevent=200 extraspread=40 spangti=no

echo "Packaging Logs"

# Pack and cleanup logs
7z a "$OBSID"_rsl_logs.7z *.log
rm *.log

echo "Cleaning up temporary files..."

# Enter the directory, remove duplicate files and gzip
rm "$OBSID"_rsl_reproc/*.hk
echo "GZipping..."
gzip "$OBSID"_rsl_reproc/*

echo "Running Xtend Pipeline..."
# Run the pipeline
xapipeline indir="$OBSID" outdir="$OBSID"_xtd_reproc steminputs=xa"$OBSID" stemoutputs=DEFAULT entry_stage=1 exit_stage=3 instrument=xtend verify_input=no create_ehkmkf=no calc_pointing=yes calc_optaxis=yes seed=650504

echo "Packaging Logs"

# Pack and cleanup logs
7z a "$OBSID"_xtd_logs.7z *.log
rm *.log

echo "Cleaning up temporary files..."

# Enter the directory, remove duplicate files and gzip
rm "$OBSID"_xtd_reproc/*.hk
echo "GZipping..."
gzip "$OBSID"_xtd_reproc/*

cd ..

echo "Setting up Analysis Directory..."

mkdir analysis

cp "$OBSID"_xtd_reproc/*_cl.evt.gz analysis/
cp "$OBSID"_rsl_reproc/*_cl.evt.gz analysis/
cp "$OBSID"_rsl_reproc/xa"$OBSID".ehk.gz "$OBSID"_rsl_reproc/xa"$OBSID"rsl_px1000_exp.gti.gz "$OBSID"_xtd_reproc/xa"$OBSID"xtd_p031100010.bimg.gz analysis/

echo "Complete!"
