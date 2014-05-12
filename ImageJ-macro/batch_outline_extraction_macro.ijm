//---- IMPORTANT: for this macro to work, the color setting must be follow:
//---- Foreground: black, background: white (the default)
//---- to set, go to Edit> Options> Colors

// Author: Wong Jin Yung

//-------- get largest particle, from KeepLargestPrticle.txt(author: G. Landini), updated march 2013--------//
function getlargest (){
  run("Set Measurements...", "area");
        run("Analyze Particles...", "minimum=0 maximum=999999999999999 bins=20 show=Nothing clear record");
        
        ar=0;
        for (j=0; j < nResults; j++) {    
            ca = getResult('Area', j);
            if (ca>ar) ar=ca;
        }
        
        for (j=0; j < nResults; j++) {
            x = getResult('XStart', j);
            y = getResult('YStart', j);
            ca = getResult('Area', j);
            if (ca<ar) {
                doWand(x,y);
                run("Fill");
                    }
        }
        run("Select None");
}
//-------- get largest particle end--------//

macro "Batch Outline Extraction"{
  // get the argument from shell, passed from runmacro() in R
        folders = getArgument;
        if (folders == "") 
          exit ("No argument!"); 
        delimiter = "*";
        parts= split(folders, delimiter); // give array of str separated by delimiter
        openDir = parts[0]; // take the first item of parts array as input dir
        saveDir = parts[1]; // and the second as output dir
        openList= getFileList(openDir); // array object containing all the file names
        
        setBatchMode(true); //images are not updated during batch mode

        for (i=0; i< lengthOf(openList); i++) { 
                open(openDir + openList[i]);
                run("8-bit");
                run("Enhance Contrast", "saturated=0.1"); // saturated value higher= higher contrast
                run("Smooth");
                savePath=saveDir + openList[i];

                // remove file extension from 'openList', so a suffix to the file name can be added.
                dotIndex = lastIndexOf(savePath, "."); 
                if (dotIndex!=-1) 
                        savePath = substring(savePath, 0, dotIndex); 
                // set threshold
                setThreshold(30, 255);
                // remove the dirts
    						getlargest();
                run("Create Mask");
                run("Fill Holes");
                 run("Duplicate...", "title=masktooutline");
                 run("Erode");
                getlargest();                
                mask = getImageID();
                run("Outline");
                run("Invert");
                FractalID = getImageID();
                run("Options...", "iterations=1 count=1 black edm=32-bit do=Nothing");
                run("Voronoi");
                run("8-bit");
                run("Invert");
    						// save the outline as image
                saveAs("Tiff", savePath + "_outline");

                // savexy coords, a pixel wide for outline
                argument= "background=255 suppress save=[" + savePath + "_xycoords.txt]";
                run("Save XY Coordinates...", argument);
                
                // get Fractal Dimension, updated March 2013
                selectImage(FractalID);
                run("Fractal Box Count...", "box=2,3,4,6,8,12,16,32,64,256,512");
                D = getResult("D");
                selectWindow("Results");
                        run("Close"); 
                selectWindow("Plot");
                             run("Close"); 

                // get dimentionless shape indices, using Particle8 plugin by Landini, G.
                // http://www.dentistry.bham.ac.uk/landinig/software/software.html
                selectImage(mask);
                run("Particles8 ", "white exclude label morphology show=Particles minimum=0 maximum=9999999 redirect=None");
                setResult("FracD", 0, D ); 
                setResult("Label", 0, openList[i]);
                saveAs("Results", savePath + "_shapedescriptors.txt");
                selectWindow("Results");
                        run("Close");
                       
                // Clean-up of any image opened by ImageJ
                while (nImages() > 0) {
                           selectImage(nImages());  
                        run("Close");
                }

                // Progress window
                if (i == 0){
                        title = "[Progress]"; // must include "[]"!
                        run("Text Window...", "name="+ title +" width=25 height=1"); // open new text window
                }
                k= i + 1;        
                print(title, "\\Update:" + k + "/" + openList.length + " images (" + round(k / openList.length*100)+"%)"); // update the progress window and show the progress
                
                }
                print(title, "\\Close"); // close the progress window
}  