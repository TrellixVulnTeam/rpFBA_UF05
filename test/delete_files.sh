
declare -a arr=("rp_14_1.sbml.xml" "rp_26_1.sbml.xml" "rp_19_1.sbml.xml" "rp_20_1.sbml.xml" "rp_18_1.sbml.xml" "rp_21_1.sbml.xml" "rp_5_2.sbml.xml" "rp_11_1.sbml.xml" "rp_9_1.sbml.xml" "rp_25_1.sbml.xml" "rp_24_2.sbml.xml" "rp_1_2.sbml.xml" "rp_4_1.sbml.xml" "rp_11_2.sbml.xml" "rp_29_1.sbml.xml" "rp_30_1.sbml.xml" "rp_8_2.sbml.xml" "rp_15_1.sbml.xml" "rp_28_1.sbml.xml" "rp_13_1.sbml.xml" "rp_23_1.sbml.xml" "rp_8_1.sbml.xml" "rp_5_1.sbml.xml" "rp_17_1.sbml.xml" "rp_24_1.sbml.xml" "rp_23_2.sbml.xml" "rp_22_1.sbml.xml" "rp_6_1.sbml.xml" "rp_7_1.sbml.xml" "rp_10_1.sbml.xml" "rp_27_1.sbml.xml" "rp_12_1.sbml.xml" "rp_16_2.sbml.xml" "rp_25_2.sbml.xml" "rp_16_1.sbml.xml")

cp test_input.tar test_input.tar.xz

arraylength=${#arr[@]}

echo $arraylength

for i in "${arr[@]}"
	do
		unxz < test_input.tar.xz | tar --delete $i | xz > test_input1.tar.xz
		mv test_input1.tar.xz test_input.tar.xz
	done

