#!/bin/bash
################################################################################
#
# PREPARATIONS
#
################################################################################

export LC_NUMERIC="en_US.UTF-8"

################################################################################
# Simple formatting

bold=$(tput bold)
normal=$(tput sgr0)

################################################################################


function Help() {
    cat <<HELP

Usage:

$(basename $0) ${bold}-s${normal} full path to SWI ${bold} ${bold}-t${normal} full path to TOF ${bold} ${bold}-u${normal} path to hippunfold output

--------------------------------------------------------------------------------
Input arguments:

    -s: full path to processed SWI nifti. It should be vessel masked and resliced to anatomical space
    -t: full path to processed TOF nifti. It should be vessel masked and resliced to anatomical space
    -u: path to hippunfold output containing anat,coords,surf, etc.

--------------------------------------------------------------------------------

Example:

$(basename $0) -s SWIpath/rSWI_vessel_masked.nii -t TOFpath/rTOF_vessel_masked.nii -u hippunfoldpath/s-7553

This command is a wrapper to wb_command -volume-to-surface-mapping. 
--------------------------------------------------------------------------------
Script was created by: Viktor Pfaffenrot (01-2024)
--------------------------------------------------------------------------------
HELP
    exit 1
}

if [[ "$1" == "-h" || $# -eq 0 ]]; then
    Help >&2
fi



################################################################################
#
# PARSE INPUT ARGUMENTS
#
################################################################################

while getopts "s:t:u:" OPT; do
    case $OPT in
    h) #help
        Help
        exit 0
        ;;
    s) # input file path
        SWI=$OPTARG
        ;;
    t) # input file path
        TOF=$OPTARG
        ;;        
    u) # input file path
        hippunfold_path=$OPTARG
        ;;        
    \?) # report error
        echo "$HELP" >&2
        exit 1
        ;;
    esac
done


################################################################################
#
# DEFAULTS
#
################################################################################
interpolation_tpye=cubic
hemis=("L" "R")
SURFS=("inner" "outer")



for hemi in ${hemis[@]}; do
	for SURF in ${SURFS[@]}; do
		wb_command -volume-to-surface-mapping  ${SWI} \
		"${hippunfold_path}/surf/"$(echo *hemi-${hemi}_space-T2w_den-0p5mm_label-hipp_${SURF}.surf.gii) \
		$(dirname "$SWI")/SWI_to_${SURF}_${hemi}.shape.gii \
		-$interpolation_tpye
		
		wb_command -volume-to-surface-mapping  ${TOF} \
		"${hippunfold_path}/surf/"$(echo *hemi-${hemi}_space-T2w_den-0p5mm_label-hipp_${SURF}.surf.gii) \
		$(dirname "$TOF")/TOF_to_${SURF}_${hemi}.shape.gii \
		-$interpolation_tpye		
	done
done



























