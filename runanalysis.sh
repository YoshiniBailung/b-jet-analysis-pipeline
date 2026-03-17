#!/bin/bash
python run.py

set -e

ANALYSIS_DIR=$(pwd)

HIST_DIR="$ANALYSIS_DIR/outputs/histograms"
FIT_DIR="$ANALYSIS_DIR/outputs/fits"

MACRO_SRC_1D="$ANALYSIS_DIR/macros/TemplateFitter.C"
MACRO_SRC_2D="$ANALYSIS_DIR/macros/TemplateFitter2D.C"
MACRO_SRC_DB="$ANALYSIS_DIR/macros/TemplateFitterDb.C"

alienv enter O2Physics/latest-master-o2 -c bash <<EOF

echo "Entered O2Physics environment"

for PTBIN in \$(ls $HIST_DIR)
do
    echo "Processing pT bin: \$PTBIN"

    BIN_HIST_DIR="$HIST_DIR/\$PTBIN"
    BIN_FIT_DIR="$FIT_DIR/\$PTBIN"

    mkdir -p "\$BIN_FIT_DIR"

    cp "\$BIN_HIST_DIR"/*.root "\$BIN_FIT_DIR/"


    cp "$MACRO_SRC_1D" "\$BIN_FIT_DIR/"
    cp "$MACRO_SRC_2D" "\$BIN_FIT_DIR/"
    cp "$MACRO_SRC_DB" "\$BIN_FIT_DIR/"

    cd "\$BIN_FIT_DIR"

	root -l -b -q "TemplateFitter.C(\"data_1d.root\",\"hist\",\"template_1d_b.root\",\"hist_b\",\"template_1d_c.root\",\"hist_c\",\"template_1d_lf.root\",\"hist_lf\",\"\$PTBIN\")"

	root -l -b -q "TemplateFitter2D.C(\"data_2d.root\",\"hist\",\"template_2d_b.root\",\"hist_b\",\"template_2d_c.root\",\"hist_c\",\"template_2d_lf.root\",\"hist_lf\",\"\$PTBIN\")"

	root -l -b -q "TemplateFitterDb.C(\"data_Db.root\",\"hist\",\"template_Db_b.root\",\"hist_b\",\"template_Db_c.root\",\"hist_c\",\"template_Db_lf.root\",\"hist_lf\",\"\$PTBIN\")"

    echo "Finished pT bin: \$PTBIN"

    cd "$ANALYSIS_DIR"

done

echo "All fits completed"

EOF
