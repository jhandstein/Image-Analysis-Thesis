print("\\Clear");

//choosing directory where input images are saved
inputdir = getDirectory("Choose directory containing Yokogawa CQ1 raw images...")
//generate array with each file name(?)
list=getFileList(inputdir); 

//create string with unique identifier of experiment. text in brackets should be changed depending on experiment
experimentID = "210121_SEM6_DRAQ5\\"
print("Processing " + experimentID + " ...")

//creating folder for output images
outputdir = getDirectory("Select  2_Data Extraction...")
outputdir = outputdir + experimentID ;
File.makeDirectory(outputdir);

//position of first image file in directory (in case of creating new analysis folder in raw image folder this has to be 1)
pos = 0

//open one specific channel for each field. choice of C4 is random. this is in purpose of counting the amount of fields. "n" is the amount of fields/number of slices in stack
run("Image Sequence...", "open=["+inputdir+list[pos]+"] file=C4 convert sort");
n = nSlices
close();

for (i = 0; i < n; i++) {

	iset = i+1; 	//number of processed image sets
	print("Image Set: "+iset+"");

    //get string for the name of current field
    open(""+inputdir+list[pos]+"");
	field = getTitle();
	//fetch unique identifier for each field
	field = substring(field, 0, 10);
	print(field);
	close();

	fielddir = outputdir + "/"+field+"/";
	File.makeDirectory(fielddir);

	run("Image Sequence...", "open=["+inputdir+list[pos]+"] file="+field+" convert sort");
	rename("stack_"+field+"");
	stack = getTitle();

	//get rid of brightfield channel. set up stack for further analysis
	setSlice(5);
	run("Delete Slice");
	setSlice(4);

	run("Duplicate...", "title=mask");
	//image=getTitle();
	run("Set Scale...", "distance=0 known=0 unit=pixel");

	close(stack);

	//run both filters with a fixed value
	run("Median...", "radius=5");
	run("Kuwahara Filter", "sampling=5");
	setOption("ScaleConversions", true);

	setAutoThreshold("MinError dark");
	//setThreshold(25, 255);
	run("Convert to Mask");
	run("Open");
	run("Close-");
	run("Watershed");

	//since objects at the edge still got detected, make rectangle with a distance of one pixel away from all edges and create new image with it
	makeRectangle(1, 1, 2558, 2158);
	run("Duplicate...", "title=mask_cropped");
	
	setting = "1_watershed";
	saveAs("Tiff", ""+fielddir+""+setting+"_"+field+".tif");
	image=getTitle();
	close("mask");
	
	//save mask with ROI outlines
	selectWindow(""+image+"");
	run("Duplicate...", "title=mask-1");
	run("Analyze Particles...", "size=150-Infinity display exclude clear add");
	roiManager("Show All without labels");

	setting = "2_ROI_outlines";
	saveAs("Tiff", ""+fielddir+""+setting+"_"+field+".tif");
	close();

	
	//create an image where every pixel inside a cell is a value of 1 and every pixel outside is a pixel of 0
	selectWindow(""+image+"");
	run("Divide...", "value=255");
	setting = "3_mask";
	saveAs("Tiff", ""+fielddir+""+setting+"_"+field+".tif");
	mask = getTitle();


		for (k = 0; k < 3; k++) {
		//opening the first the fluorescent protein channel (Channel i+1). count is position of the respective channel. count_t is array with number of the three quantative channels
		count = pos + k;
		count_t = newArray(1, 2, 3);
		open(""+inputdir+list[count]+"");
		name_uncr=getTitle();

		//make image with same dimensions as cropped mask
		makeRectangle(1, 1, 2558, 2158);
		run("Duplicate...", "title="+name_uncr+"_cr");
		rename("Ch_"+count_t[k]+"");
		name_Ch=getTitle();
		run("Set Scale...", "distance=0 known=0 unit=pixel");
		//close original image
		close(name_uncr);

		imageCalculator("Multiply create", ""+name_Ch+"",""+mask+"");
		another_byproduct=getTitle();

		//make image with same dimensions as cropped mask
		makeRectangle(1, 1, 2556, 2156);
		run("Duplicate...", "title=doublecropped");
		close(another_byproduct);
		
		//setAutoThreshold("MinError dark");
		setThreshold(1, 65536);
		run("Analyze Particles...", "size=150-Infinity display exclude clear add");
		saveAs("Results", ""+outputdir+""+field+"_Ch"+count_t[k]+"_intensities.csv");
		close("Results");

		setting = "objects_"+count_t[k]+"";
		saveAs("Tiff", ""+fielddir+""+setting+".tif");

		close();
		close(name_Ch);
		}
	
	close(mask);
    
    pos = pos + 5; 
}
