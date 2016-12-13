package glcm;

// glcm.GLCMTexture.java v0.011
// Writed by Angel Zeitoune, based on previews versions
// 06/04/2016
// https://github.com/Crandelz/Texture-Analysis/
//
// Kota Miura (miura@embl.de), Toby C. Cornish, Julio E. Cabrera
// 05/24/2011
//=================================================================================================
//  	Modified from GLCM_TextToo.java v0.008 by Toby C. Cornish for use of the class as library.
//	 	See documentation by Toby copied from GLCM_TextToo.java.
//
//		Following changes were made:
//			- actual processing and math parts were moved from GLCM_TextToo.java to this code.
//			-
//	 		- added several constructors
//	 		- calcGLCM(ImagePlus imp) calculate GLCM
//	 		- added access to Texture parameters by getResultsArray()
//			- offsetX and offsetY calculated by trigonometry
//	 		- added "reset Results table" check box in the dialog, in method showDialog()
//	 		- for use from scripting languages, see example below.
//
//		example javascript
//			importClass(Packages.emblcmci.glcm.GLCMtexture);
//			g = new GLCMtexture(IJ.getImage(), 2, 45, true, true);
//			glcm = g.calcGLCM();
//			ht = g.getResultsArray();
//			IJ.log(ht.get("Contrast"));
//
//=================================================================================================
//GLCM_Texture_Too, v. 0.008
//Toby C. Cornish, tcornis3@jhmi.edu (or toby@tobycornish.com)
//11/26/2007
//modified from GLCM_Texture (Gray Level Correlation Matrix Texture Analyzer) v0.4, Author: Julio E. Cabrera, 06/10/05
//
//CHANGELOG:
//------------------------------------------------------------------------------------------------
//11/26/07 GLCM_TextureToo v0.008 - minor fixes; added references section -- tc
//11/22/07 GLCM_TextureToo v0.007 - minor fixes to parameters -- tc
//                              - moved mean and stdev calculations to common area -- tc
//11/20/07 GLCM_TextureToo v0.006 - confirmed results against parker and Xite -- tc
//                              - added preliminary shade, prominence, variance, homogeneity, inertia
//                              - differentiated Haralick correlation from Walker correlation
//11/19/07 GLCM_TextureToo v0.005 - changed method of calculating correlation -- tc
//                              - should be closest in nomenclature to Walker, et al.  1995
//11/19/07 GLCM_TextureToo v0.004 - corrected column name of idm -- tc
//11/13/07 GLCM_TextureToo v0.003 - changed from roi.contains() to byte[] checking -> much faster -- tc
//11/13/07 GLCM_TextureToo v0.002 - added progress bar --tc
//11/11/07 GLCM_TextureToo v0.001 - inherited portions of the codebase from GLCM_Texture 0.4 -- tc
//                              - fundamental rewrite of GLCM calculation
//                              - added irregular ROI support
//                              - added symmetrical/non-symmetrical GLCM calculations
//                              - corrected directionality (phi) to be 0,45,90,135
//
//=================================================================================================
//
//References:
//R.M. Haralick, Texture feature for image classification, IEEE Trans. SMC 3 (1973) (1), pp. 610-621.
//Conners, R.W., Trivedi, M.M., and Harlow, C.A., Segmentation of a High-Resolution Urban Scene
//  Using Texture Operators, CVGIP(25), No. 3, March, 1984, pp. 273-310.
//Walker, RF, Jackway, P and Longstaff, ID (1995) Improving Co-occurrence Matrix Feature Discrimination.'
//  In  DICTA '95, 3rd Conference on Digital Image Computing: Techniques and Application, 6 - 8 December,
//  1995, pages 643-648.
//Parker, JR, Algorithms for Image Processing and Computer Vision, John Wiley & Sons, 1997.
//Image processing lab, Department of Informatics, University of Oslo. Xite v1.35: glcmParameter.c, v1.30
//  2004/05/05 07:34:19 (2004)
//
//Bankman, IN, Ed., Handbook of Medical Image Processing and Analysis, 2nd Edition, 2009. Academic Press


import java.awt.Rectangle;

import ij.IJ;
import ij.process.ImageProcessor;

public class GLCMTexture {

    // Parameters

    public int pixelOffset = 1;             // step (displacement or distance) between reference and the neighbour pixel
    public int phi = 0;                     // Angle
    public boolean symmetry = true;         // Perform symmetry calculation
    public boolean doAngularSecondMoment = true;
    public boolean doContrast = true;
    public boolean doCorrelation = true;
    public boolean doInverseDifferenceMoment = true;
    public boolean doEntropy = true;
    public boolean doEnergy = true;
    public boolean doInertia = true;
    public boolean doHomogeneity = true;
    public boolean doProminence = true;
    public boolean doVariance = true;
    public boolean doShade = true;

    // Internal common statistical parameters
    private double [][] glcm;
    private double meanx=0.0;
    private double meany=0.0;
    private double stdevx=0.0;  // really is the variance
    private double stdevy=0.0;

    /*private double [] meanx1 = new double [256];
    private double [] meany1 = new double [256];
    private double [] varianceX = new double [256];
    private double [] varianceY = new double [256];*/



    public double AngularSecondMoment;
    public double InverseDifferenceMoment;
    public double Contrast;
    public double Energy;
    public double Entropy;
    public double Homogeneity;
    public double Variance;
    public double Shade;
    public double Prominence;
    public double Inertia;
    public double Correlation;
    public double GLCMsum;

    public String[] paramA =
            {
                    "Angular Second Moment",
                    "Inverse Difference Moment",
                    "Contrast",
                    "Energy",
                    "Entropy",
                    "Homogeneity",
                    "Variance",
                    "Shade",
                    "Prominence",
                    "Inertia",
                    "Correlation",
                    "Sum of all GLCM elements"
            };

    /**Constructor for use as library
     * all the parameters will be by default be measured (true).
     * @param d
     * @param phi
     * @param symmetry
     */

    public GLCMTexture(){
    }

    private void doBasicStats(){
        double [] px = new double [256];
        double [] py = new double [256];
        meanx=0.0;
        meany=0.0;
        stdevx=0.0;
        stdevy=0.0;

        // Px(i) and Py(j) are the marginal-probability matrix; sum rows (px) or columns (py)
        // First, initialize the arrays to 0
        for (int i=0;  i<256; i++){
            px[i] = 0.0;
            py[i] = 0.0;

            /*meanx1[i] = 0.0;
            meany1[i] = 0.0;
            varianceX[i] = 0.0;
            varianceY[i] = 0.0;*/
        }

        // sum the glcm rows to Px(i)
        for (int i=0;  i<256; i++) {
            for (int j=0; j<256; j++) {
                px[i] += glcm [i][j];
            }
            //meanx1[i] = i*px[i];
            //IJ.log("meanx1 : " + meanx1[i]);
        }

        // sum the glcm rows to Py(j)
        for (int j=0;  j<256; j++) {
            for (int i=0; i<256; i++) {
                py[j] += glcm [i][j];
            }
            //meany1[j] = j*py[j];
            //IJ.log("meany1 : " + meany1[j]);
        }

        // calculate meanx and meany
        for (int i=0;  i<256; i++) {

            meanx += (i*px[i]);
            meany += (i*py[i]);
        }

        // calculate stdevx and stdevy
        for (int i=0;  i<256; i++) {
            stdevx += ((i-meanx)*(i-meanx)*px[i]);
            stdevy += ((i-meany)*(i-meany)*py[i]);

            //varianceX[i] += ((i-meanx1[i])*(i-meanx1[i])*py[i]);
            //varianceY[i] += ((i-meany1[i])*(i-meany1[i])*px[i]);

            //IJ.log("varianceX: " + varianceX[i]);
            //IJ.log("varianceY: " + varianceY[i]);
        }
    }

    //=====================================================================================================
    // This is the generic moments function from parker -- may implement in the future
    // k is the power for the moment

	/*
if (doMoments == true){
	double y=0.0;
	double z;
	double k;
	double moments = 0.0;
	for (int i=0;  i<256; i++)  {
		for (int j=0; j<256; j++) {
			if (k>0) {
				z = Math.pow ((i-j), k);
			} else {
				if (i == j) continue;
				z = Math.pow ((i-j), -1*k);
				z = glcm[i][j]/z;
			}
			moments += z * glcm[i][j];
		}
	}
	rt.setValue("Angular Second Moment", row, asm);
}
	 */
    //=====================================================================================================
    // calculate the angular second moment (asm)
    // also known as 'energy' (formula 15.38, Bankman, 2009)

    public double getAngularSecondMoment(){
        double asm = 0.0;
        for (int i=0;  i<256; i++)  {
            for (int j=0; j<256; j++) {
                asm += (glcm[i][j]*glcm[i][j]);
            }
        }
        return asm;
    }
    //===============================================================================================
    // calculate the inverse difference moment (idm) (Walker, et al. 1995)
    // this is calculated using the same formula as Conners, et al., 1984 "Local Homogeneity"
    // (formula 15.40, Bankman, 2009)

    public double getInverseDifferenceMoment(){
        double IDM = 0.0;
        for (int i=0;  i<256; i++)  {
            for (int j=0; j<256; j++) {
                IDM += ((1/(1+(Math.pow(i-j,2))))*glcm[i][j]);
            }
        }
        return IDM;
    }

    //===============================================================================================
    // (formula 15.39, Bankman, 2009) energy weighted by pixel value difference
    public double getContrast(){
        double contrast=0.0;

        for (int i=0;  i<256; i++)  {
            for (int j=0; j<256; j++) {
                contrast += ((i-j)*(i-j) * (glcm[i][j]));
            }
        }
        return contrast;
    }
    //===============================================================================================
    // calculate the energy
    // - same as Angular 2nd Moment (see above), so may be delete this.
    // - TODO (in calgary website, ASM is square rooted to give energy)
    public double getEnergy(){
        double energy = 0.0;
        for (int i=0;  i<256; i++)  {
            for (int j=0; j<256; j++) {
                energy += Math.pow(glcm[i][j],2);
            }
        }
        return energy;
    }
    //===============================================================================================
    // calculate the entropy (Haralick et al., 1973; Walker, et al., 1995)
    // -TODO there are some difference in textbooks as well.
    public double getEntropy(){
        double entropy = 0.0;
        for (int i=0;  i<256; i++)  {
            for (int j=0; j<256; j++) {
                if (glcm[i][j] != 0) {
                    entropy = entropy-(glcm[i][j]*(Math.log(glcm[i][j])));
                    //the next line is how Xite calculates it -- I am not sure why they use this, I do not think it is correct
                    //(they also use log base 10, which I need to implement)
                    //entropy = entropy-(glcm[i][j]*((Math.log(glcm[i][j]))/Math.log(2.0)) );
                }
            }
        }
        return entropy;
    }
    //===============================================================================================
    // calculate the homogeneity (Parker)
    // "Local Homogeneity" from Conners, et al., 1984 is calculated the same as IDM above
    // Parker's implementation is below; absolute value of i-j is taken rather than square
    // - matlab textbook also uses non-squred absolute difference |i-j|
    // -- using absolute value, flat image (diagonal) will be 1.
    public double getHomogeneity(){
        double homogeneity = 0.0;
        for (int i=0;  i<256; i++) {
            for (int j=0; j<256; j++) {
                homogeneity += glcm[i][j]/(1.0+Math.abs(i-j));
            }
        }
        return homogeneity;
    }

    //===============================================================================================
    // calculate the variance ("variance" in Walker 1995; "Sum of Squares: Variance" in Haralick 1973)

    public double getVariance(){
        double variance = 0.0;
        double mean = 0.0;

        mean = (meanx + meany)/2;
		/*
	// this is based on xite, and is much greater than the actual mean -- it is here for reference only
	for (int i=0;  i<256; i++)  {
		for (int j=0; j<256; j++) {
			mean += glcm[i][j]*i*j;
		}
	}
		 */

        for (int i=0;  i<256; i++)  {
            for (int j=0; j<256; j++) {
                variance += (Math.pow((i-mean),2)* glcm[i][j]);
            }
        }
        return variance;
    }

    /** Shade
     * calculate the shade (Walker, et al., 1995; Connors, et al. 1984)
     * @return
     */
    public double getShade(){
        double shade = 0.0;

        // calculate the shade parameter
        for (int i=0;  i<256; i++) {
            for (int j=0; j<256; j++) {
                shade += (Math.pow((i+j-meanx-meany),3)*glcm[i][j]);
            }
        }
        return shade;
    }

    //==============================================================================================
    // calculate the prominence (Walker, et al., 1995; Connors, et al. 1984)
    public double getProminence(){
        double prominence=0.0;

        for (int i=0;  i<256; i++) {
            for (int j=0; j<256; j++) {
                prominence += (Math.pow((i+j-meanx-meany),4)*glcm[i][j]);
            }
        }
        return prominence;
    }

    //===============================================================================================
    // calculate the inertia (Walker, et al., 1995; Connors, et al. 1984)
    public double getInertia(){
        double inertia = 0.0;
        for (int i=0;  i<256; i++)  {
            for (int j=0; j<256; j++) {
                if (glcm[i][j] != 0) {
                    inertia += (Math.pow((i-j),2)*glcm[i][j]);
                }
            }
        }
        return inertia;
    }
    //=====================================================================================================
    /** calculate the correlation
     *  methods based on Haralick 1973 (and MatLab), Walker 1995 are included below
     * Haralick/Matlab result reported for correlation currently; will give Walker as an option in the future
     */
    public double getCorrelation(){
        double correlation=0.0;
        double correlation2=0.0;

        // calculate the correlation parameter
        for (int i=0;  i<256; i++) {
            for (int j=0; j<256; j++) {
                correlation += (((i - meanx) * (j - meany)) * glcm[i][j]);

                /*if (varianceX[i] != 0.0 && varianceY[j] != 0.0) {
                    correlation2 += (((i - meanx1[i]) * (j - meany1[j])) * glcm[i][j]) / Math.sqrt(varianceX[i] * varianceY[j]);
                }*/
            }
        }
        correlation /= Math.sqrt(stdevx * stdevy);

        //IJ.log("correlation 1 : " + correlation);
        //IJ.log("correlation 2 : " + correlation2);
        return correlation;
    }

    //===============================================================================================
    // calculate the sum of all glcm elements
    public double getGLCMsum(){
        double sum = 0.0;
        for (int i=0; i<256; i++)  {
            for (int j=0; j<256; j++) {
                sum = sum + glcm[i][j];
            }
        }
        return sum;
    }
    //=========================================================================================

    public void CalculateStatistics(){
        doBasicStats();

        AngularSecondMoment = getAngularSecondMoment();
        InverseDifferenceMoment = getInverseDifferenceMoment();
        Contrast = getContrast();
        Energy = getEnergy();
        Entropy = getEntropy();
        Homogeneity = getHomogeneity();
        Variance = getVariance();
        Shade = getShade();
        Prominence = getProminence();
        Inertia = getInertia();
        Correlation = getCorrelation();
        GLCMsum = getGLCMsum();
    }

    /**main part that does the calculation of GLCM
     *
     * @param ip
     */
    public double [][] calcGLCM(ImageProcessor ip){
        // use the bounding rectangle ROI to roughly limit processing
        Rectangle roi = ip.getRoi();

        // get byte arrays for the image pixels and mask pixels
        int width = ip.getWidth();      // window width

        byte [] pixels = (byte []) ip.getPixels();
        byte [] mask = ip.getMaskArray();

        //double totalPixels = roi.height * roi.width;
        //if (symmetry) totalPixels = totalPixels * 2;
        //double pixelProgress = 0;
        double pixelCount = 0;

        //IJ.log("roi width: " + roi.width);
        //IJ.log("roi height: " + roi.height);
        //IJ.log("ip width: " + ip.getWidth());
        //IJ.log("ip height: " + ip.getHeight());

        //====================================================================================================
        // compute the Gray Level Correlation Matrix

        //int offsetX = 1;
        //int offsetY = 0;

        // set our offsets based on the selected angle
        //		if (phi == 0) {
        //			offsetX = d;
        //			offsetY = 0;
        //		} else if (phi == 45) {
        //			offsetX = d;
        //			offsetY = -d;
        //		} else if (phi == 90) {
        //			offsetX = 0;
        //			offsetY = -d;
        //		} else if (phi == 135) {
        //			offsetX = -d;
        //			offsetY = -d;
        //		} else {
        //			// the angle is not one of the options
        //			IJ.showMessage("The requested angle,"+phi+", is not one of the supported angles (0,45,90,135)");
        //		}
        double rad = Math.toRadians(-1.0 * phi);
        int offsetX = (int) (pixelOffset * Math.round(Math.cos(rad)));
        int offsetY = (int) (pixelOffset * Math.round(Math.sin(rad)));
        //IJ.log("Angle: " + phi + " - offset X: " + offsetX + " - offset Y: " + offsetY );

        int value;      // value at pixel of interest
        int dValue;     // value of pixel at offset
        glcm = new double [256][256];   // GLCM matrix. Only for 8-bit images

        // loop through the pixels in the ROI bounding rectangle
        for (int y=roi.y; y<(roi.y + roi.height); y++) 	{
            for (int x=roi.x; x<(roi.x + roi.width); x++)	 {
                // check to see if the pixel is in the mask (if it exists)
                if ((mask == null) || ((0xff & mask[(((y-roi.y)*roi.width)+(x-roi.x))]) > 0) ) {
                    // check to see if the offset pixel is in the roi
                    int dx = x + offsetX;
                    int dy = y + offsetY;
                    if ( ((dx >= roi.x) && (dx < (roi.x+roi.width))) && ((dy >= roi.y) && (dy < (roi.y+roi.height))) ) {
                        // check to see if the offset pixel is in the mask (if it exists)
                        if ((mask == null) || ((0xff & mask[(((dy - roi.y) * roi.width) + (dx - roi.x))]) > 0)) {
                            value = 0xff & pixels[(y * width) + x];
                            dValue = 0xff & pixels[(dy * width) + dx];
                            glcm[value][dValue]++;
                            pixelCount++;
                        }
                    }

                    // if symmetry is selected, invert the offsets and go through the process again
                    if (symmetry) {
                        dx = x - offsetX;
                        dy = y - offsetY;
                        if ( ((dx >= roi.x) && (dx < (roi.x+roi.width))) && ((dy >= roi.y) && (dy < (roi.y+roi.height))) ) {
                            // check to see if the offset pixel is in the mask (if it exists)
                            if ((mask == null) || ((0xff & mask[(((dy-roi.y)*roi.width)+(dx-roi.x))]) > 0) ) {
                                value = 0xff & pixels[(y*width)+x];
                                dValue = 0xff & pixels[(dy*width) + dx];
                                glcm [dValue][value]++;
                                pixelCount++;
                            }
                        }
                    }
                }
                //pixelProgress++;
                //IJ.showProgress(pixelProgress/totalPixels);
            }
        }

        // convert the GLCM from absolute counts to probabilities
        for (int i=0; i<256; i++)  {
            for (int j=0; j<256; j++) {
                glcm[i][j] = (glcm[i][j])/(pixelCount);
            }
        }
        return glcm;
    }

}