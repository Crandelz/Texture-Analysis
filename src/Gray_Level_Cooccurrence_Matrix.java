/**
 * Created by Angel on 06/04/2016.
 */

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.measure.ResultsTable;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import glcm.GLCMTexture;

public class Gray_Level_Cooccurrence_Matrix implements PlugInFilter {

    boolean processCurrentChannel = true;
    public boolean newInstance = true;

    private GLCMTexture gl = new GLCMTexture();


    public int setup(String arg, ImagePlus imp) {
        //IJ.log("new instance " + newInstance);
        if (imp!=null && !showDialog()) return DONE;
        newInstance = false;
        return DOES_8G+DOES_STACKS+SUPPORTS_MASKING;
    }

    public void run(ImageProcessor ip) {
        if (!processCurrentChannel || ip.getSliceNumber() == IJ.getImage().getChannel()) {
            // Only process selected channel
            gl.calcGLCM(ip);
            gl.CalculateStatistics();
            writeToResultsTable(true);
        }
    }

    public void writeToResultsTable(boolean ShowResults){

        ResultsTable rt = ResultsTable.getResultsTable();
        int row = rt.getCounter();
        rt.incrementCounter();
        if (gl.doAngularSecondMoment)
            rt.setValue("Angular Second Moment", row, gl.AngularSecondMoment);
        if (gl.doInverseDifferenceMoment)
            rt.setValue("Inverse Difference Moment", row, gl.InverseDifferenceMoment);
        if (gl.doContrast)
            rt.setValue("Contrast", row, gl.Contrast);
        if (gl.doEnergy)
            rt.setValue("Energy", row, gl.Energy);
        if (gl.doEntropy)
            rt.setValue("Entropy", row, gl.Entropy);
        if (gl.doHomogeneity)
            rt.setValue("Homogeneity", row, gl.Homogeneity);
        if (gl.doVariance)
            rt.setValue("Variance", row, gl.Variance);
        if (gl.doShade)
            rt.setValue("Shade", row, gl.Shade);
        if (gl.doProminence)
            rt.setValue("Prominence", row, gl.Prominence);
        if (gl.doInertia)
            rt.setValue("Inertia", row, gl.Inertia);
        if (gl.doCorrelation)
            rt.setValue("Correlation", row, gl.Correlation);

        rt.setValue("Sum of all GLCM elements", row, gl.GLCMsum);
        if (ShowResults) rt.show("Results");
    }

    // implementation of the dialog
    private boolean showDialog() {
        GenericDialog gd = new GenericDialog("GLCM v0.011");
        String [] channels = {"Current", "All"};
        gd.addChoice("Select the channel to process", channels, "Current");
        gd.addNumericField ("Enter the size of the step in pixels",  gl.pixelOffset, 0);
        String [] angles={"0", "45", "90", "135"};
        gd.addChoice("Select the direction of the step", angles, Integer.toString(gl.phi));
        gd.addCheckbox("Symmetrical GLCM?", gl.symmetry);

        gd.addMessage("Calculate which parameters?");
        gd.addCheckbox("Angular Second Moment  ", gl.doAngularSecondMoment);
        gd.addCheckbox("Contrast  ", gl.doContrast);
        gd.addCheckbox ("Correlation  ", gl.doCorrelation);
        gd.addCheckbox ("Inverse Difference Moment  ", gl.doInverseDifferenceMoment);
        gd.addCheckbox ("Entropy   ", gl.doEntropy);
        gd.addCheckbox ("Energy   ", gl.doEnergy);
        gd.addCheckbox ("Inertia   ", gl.doInertia);
        gd.addCheckbox ("Homogeneity   ", gl.doHomogeneity);
        gd.addCheckbox ("Prominence   ", gl.doProminence);
        gd.addCheckbox ("Variance   ", gl.doVariance);
        gd.addCheckbox ("Shade   ", gl.doShade);

        gd.showDialog();
        if (gd.wasCanceled()) return false;

        processCurrentChannel = gd.getNextChoice() == "Current";
        gl.pixelOffset = (int) gd.getNextNumber();
        gl.phi = Integer.parseInt(gd.getNextChoice());
        gl.symmetry = gd.getNextBoolean();
        gl.doAngularSecondMoment = gd.getNextBoolean();
        gl.doContrast = gd.getNextBoolean();
        gl.doCorrelation = gd.getNextBoolean();
        gl.doInverseDifferenceMoment = gd.getNextBoolean();
        gl.doEntropy = gd.getNextBoolean();
        gl.doEnergy = gd.getNextBoolean();
        gl.doInertia = gd.getNextBoolean();
        gl.doHomogeneity = gd.getNextBoolean();
        gl.doProminence = gd.getNextBoolean();
        gl.doVariance = gd.getNextBoolean();
        gl.doShade = gd.getNextBoolean();

        return true;
    }
}
