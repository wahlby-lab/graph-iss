#!/usr/bin/env anduril
//$PRE startdate=$( date )
//$POST echo -e "Running time:\n$startdate -> $( date )"
//$OPT  --wrapper anduril-wrapper-slurm --threads 20
//$OPT --java-heap 4096
//$OPT -d ./result_ISS_Anduril_Pipeline_Example

import anduril.builtin._
import anduril.tools._
import anduril.sequencing._
import org.anduril.runtime._
import anduril.microarray._

object pipeline {
        val prepro_pipeline_path = INPUT(path="<GRAPH-ISS-FOLDER>/prePro_pipeline")
	prepro_pipeline_path._execute = "once"
        val predProb_pipeline_path = INPUT(path="<GRAPH-ISS-FOLDER>/prePro_pipeline/network")
	predProb_pipeline_path._execute = "once"
        val pgm_pipeline_path = INPUT(path="<GRAPH-ISS-FOLDER>/pgm_pipeline")
	pgm_pipeline_path._execute = "once"	
        val dataset_folder = INPUT(path="<DATA-FOLDER>/raw_data")
	dataset_folder._execute = "once"
        val tagList = INPUT(path="<DATA-FOLDER>/raw_data/tagList.csv")

	val img_array = Folder2Array(folder1=dataset_folder, keyMode="number", filePattern="(.*).tif", excludePattern=".*Transformed")
	val img_CSV = Array2CSV(img_array)
	val mode_percent = BashEvaluate(var1=prepro_pipeline_path, var2=img_CSV, script="""source activate pgm_pipeline
				python -W ignore -u @var1@/norm.py @var2@ 99 $( getmetadata "custom_cpu" ) @out1@ 5000
				source deactivate""")
	mode_percent._custom("cpu") ="2"

        // REGISTRATION
        val registration = BashEvaluate(var1=prepro_pipeline_path, var2=img_CSV, param1="3", param2="rigid", param3="Nuclei", param4="DO", script="""source activate pgm_pipeline
                              python -u @var1@/mainReg.py @var1@ @folder1@ @var2@ @param1@ @param2@ @param3@ @param4@ "1" "3" "noBspline"
                              source deactivate """)
        val regImg_array = Folder2Array(registration.folder1)
        val reg_CSV = Array2CSV(regImg_array)
	
	// TILING
	var tile_size_x="1330"
        var tile_size_y="980"
	val tileMap = NamedMap[BashEvaluate]("tiles")
	val tileArrayMap = NamedMap[Folder2Array]("tilesArray")
	val tileCSVMap = NamedMap[Array2CSV]("tilesCSV")
	for ((imgID, filename) <- iterArray(regImg_array.out)) {
		tileMap(imgID) = BashEvaluate(var1=regImg_array.out(imgID), param1=tile_size_x, param2=tile_size_y, script="""bfconvert -tilex @param1@ -tiley @param2@ @var1@ @folder1@/%x_%y.tif""" )
		tileMap(imgID)._custom("cpu")="4"
		tileArrayMap(imgID) = Folder2Array(folder1=tileMap(imgID).folder1, keyMode="filename")
		tileCSVMap(imgID)=Array2CSV(in=tileArrayMap(imgID))
	}
	val allTiles = makeArray(tileCSVMap)
	val allTilesJoin = CSVJoin(in=allTiles, useKeys=false)	
	val barcodesArray = NamedMap[TextFile]("barcodesArray")
	val barcodesArray_2pass = NamedMap[TextFile]("barcodesArray_2pass")
	for ((k,f)<-iterArray(allTiles).take(1)) { //take tiles idecies
		val tilesArray = CSV2Array(allTiles(k), keys="column")
		val main_files = iterArray(tilesArray)
		for ((r,file)<-main_files){
			withName(r){
				val tile = CSVFilter(in=allTilesJoin,regexp="Key=%s".format(r))
				
				// PREPROCESSING
				info("Processing tile %s.....".format(r))
				val prePro = BashEvaluate(var1=prepro_pipeline_path, var2=tile, var3=mode_percent.out1, param1="DO", script="""source activate pgm_pipeline
						python -u @var1@/main.py @var1@ @var2@ @folder1@ @out1@ @out2@ $( getmetadata "custom_cpu" ) 0.05 @var3@ @out3@ @folder2@/candidates_max.h5 @param1@
						source deactivate """) 
				prePro._custom("cpu") ="4"

				// CANDIDATE PREDICTIONS
				val predProb = BashEvaluate(var1=predProb_pipeline_path, var2=prePro.out1, script="""source activate pgm_pipeline
						python @var1@/get_proba_DNN.py @var1@ @var2@ @out1@
						source deactivate""")

				// GRAPH DECODING
				val GM = BashEvaluate(var1=pgm_pipeline_path, var2=predProb.out1, var3=tagList, var4=prePro.out2, var5=prePro.out3, var6=predProb_pipeline_path, script="""source activate pgm_pipeline
						python -u @var1@/main.py @var1@ @var2@ @out2@ @out1@ $( getmetadata "custom_cpu" ) 3 4 @out3@ blind @var3@ @var4@ @var5@ @var6@
						source deactivate""")
				GM._custom("cpu") = "1"

                                val addTileOffset = BashEvaluate(var1=GM.out1, var2=StringInput("%s".format(r)), param1=tile_size_x, param2=tile_size_y, script="""
                                                        declare $( cat @var2@ | awk '{split($0,a,"_"); print "x="a[1];print "y="a[2]}' )
                                                        echo "letters,global_X_pos,global_Y_pos,tile_ID,max_dist,seq_quality_min" > @out1@    
                                                        awk -F',' -v OFS=',' -v x=$x -v y=$y -v tile_size_x=@param1@ -v tile_size_y=@param2@ '{if (NR!=1){print $1,$2+x*tile_size_x,$3+y*tile_size_y,$4,$5,$6}}' @var1@| sed 's/"//g' >> @out1@""")
				barcodesArray(r) = addTileOffset.out1
			}
		}
	}
	val barcodes = QuickBash(in=JCSVJoin(array=barcodesArray, intersection=false), command="""sed 's/"//g' $in > $out""")
        barcodes._filename("out", "barcodes.csv")
}
