#!/usr/bin/env anduril
//$PRE startdate=$( date )
//$POST echo -e "Running time:\n$startdate -> $( date )"
//$OPT  --wrapper anduril-wrapper-slurm --threads 10
//$OPT --java-heap 4096
//$OPT -d /home/gabriele/TissueMaps_ExeDir/result_1442hippo

import anduril.builtin._
import anduril.tools._
import anduril.sequencing._
import org.anduril.runtime._
import anduril.microarray._

object pipeline {
	val prepro_pipeline_path = INPUT(path="/home/gabriele/sourceCode_repo/graph-iss/prePro_pipeline")
	prepro_pipeline_path._execute = "once"
	val predProb_pipeline_path = INPUT(path="/home/gabriele/sourceCode_repo/graph-iss/prePro_pipeline/network")
	predProb_pipeline_path._execute = "once"
	val pgm_pipeline_path = INPUT(path="/home/gabriele/sourceCode_repo/graph-iss/pgm_pipeline")
	pgm_pipeline_path._execute = "once"	
	val postpro_pipeline_path = INPUT(path="/home/gabriele/sourceCode_repo/postPro_pipeline")
	postpro_pipeline_path._execute = "once"
	val dataset_folder = INPUT(path="/data/project_data/1442_hippo")
	//dataset_folder._execute = "once"
	val tagList = INPUT(path="/home/gabriele/sourceCode_repo/graph-iss/data/taglist1441_1442.csv")

	val img_array = Folder2Array(folder1=dataset_folder, keyMode="number", filePattern="(.*).tif", excludePattern=".*Transformed")
	val img_CSV = Array2CSV(img_array)
	val mode_percent = BashEvaluate(var1=prepro_pipeline_path, var2=img_CSV, script="""
	eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
        conda activate pgm_pipeline
	python -W ignore -u @var1@/norm.py @var2@ 99 $( getmetadata "custom_cpu" ) @out1@ 50000
	conda deactivate""")
	mode_percent._custom("cpu") ="2"
        mode_percent._custom("partition") = "-p p16"

        // REGISTRATION
        val registration = BashEvaluate(var1=prepro_pipeline_path, var2=img_CSV, param1="5", param2="rigid", param3="Nuclei", param4="DO", script="""
	eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
        conda activate pgm_pipeline
        python -u @var1@/mainReg.py @var1@ @folder1@ @var2@ @param1@ @param2@ @param3@ @param4@ "512" "8" "noBspline"
        conda deactivate """)
        registration._custom("partition") = "-p p32"
        val regImg_array = Folder2Array(registration.folder1)
        val reg_CSV = Array2CSV(regImg_array)
	
	// TILING
	var tile_size_x="1024"
        var tile_size_y="1024"
	// bfconvert -tilex 512 -tiley 512 ch00_s1.ome.tiff output_tile_%x_%y.tiff	
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
		//main_files.retain((k1,v) => "%s".format(k1) == "1_1_tif") /////
		for ((r,file)<-main_files){
			withName(r){
				val tile = CSVFilter(in=allTilesJoin,regexp="Key=%s".format(r))
				
				// PREPROCESSING
				info("Processing tile %s.....".format(r))
				val prePro = BashEvaluate(var1=prepro_pipeline_path, var2=tile, var3=mode_percent.out1, param1="DO", script="""
				eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
                                conda activate pgm_pipeline
				python -u @var1@/main.py @var1@ @var2@ @folder1@ @out1@ @out2@ $( getmetadata "custom_cpu" ) 0.05 @var3@ @out3@ @folder2@/candidates_max.h5 @param1@
				conda deactivate """) 
				prePro._custom("cpu") ="4"
                                //prePro._custom("partition") = "-p p8,p16"
				//prePro._execute = "once"

				// POSTPROCESSING
				val predProb = BashEvaluate(var1=predProb_pipeline_path, var2=prePro.out1, script="""
				eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
                                conda activate pgm_pipeline
				python @var1@/get_proba_DNN.py @var1@ @var2@ @out1@
				conda deactivate""")
				//predProb._priority = 1
				//predProb._execute = "once"

				// POSTPROCESSING
				val GM = BashEvaluate(var1=pgm_pipeline_path, var2=predProb.out1, var3=tagList, var4=prePro.out2, var5=prePro.out3, var6=predProb_pipeline_path, script="""
				eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
                                conda activate pgm_pipeline
				python -u @var1@/main.py @var1@ @var2@ @out2@ @out1@ $( getmetadata "custom_cpu" ) 3 4 @out3@ blind @var3@ @var4@ @var5@ @var6@
				conda deactivate""")
				//GM._execute = "once"
				GM._custom("cpu") = "1"
                                //GM._custom("partition") = "-p p16,p8"
				//GM._priority = 2

                                val addTileOffset = BashEvaluate(var1=GM.out1, var2=StringInput("%s".format(r)), param1=tile_size_x, param2=tile_size_y, script="""
                                                        declare $( cat @var2@ | awk '{split($0,a,"_"); print "x="a[1];print "y="a[2]}' )
                                                        echo "letters,global_X_pos,global_Y_pos,tile_ID,max_dist,seq_quality_min" > @out1@    
                                                        awk -F',' -v OFS=',' -v x=$x -v y=$y -v tile_size_x=@param1@ -v tile_size_y=@param2@ '{if (NR!=1){print $1,$2+x*tile_size_x,$3+y*tile_size_y,$4,$5,$6}}' @var1@| sed 's/"//g' >> @out1@""")
				barcodesArray(r) = addTileOffset.out1
			}
		}
	}
	val barcodes = QuickBash(in=JCSVJoin(array=barcodesArray, intersection=false), command="""sed 's/"//g' $in > $out""")
}
