package mmc;

/**
 * This class provides routines for calculating roots by the combination of the Newton-Raphson method and bisection.
 * Include the function to be integrated in a class that implements the interface FunctionOfx (defined below).
 * Methods contained here are based on the Numerical Recipes (W. H. Press et al.).
 * <pre>
 * interface DFunctionOfx extends FunctionOfx{
 *     double dFunction(double x);
 * }
 * </pre>
 * For the definition of interface FunctionOfx see the manual page for class Integral.
 * @author Dmitry Chirkin
 */

public class FindRoot extends MathModel{

    private double precision;
    private int maxSteps;
    private DFunctionOfx function2use;

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes class - this constructor uses default settings
     */

    public FindRoot(){
	this(20, 1.e-6);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes class - this is the main constructor
     */

    public FindRoot(int maxSteps, double precision){
	if(maxSteps<=0){
	    Output.err.println("Warning (in Integral/Integral/1): maxSteps = "+maxSteps+" must be > 0, setting to 1");
	    maxSteps=1;
	}
	if(precision<=0){
	    Output.err.println("Warning (in Integral/Integral/2): precision = "+precision+" must be > 0, setting to 1.e-6");
	    precision=1.e-6;
	}
	this.maxSteps=maxSteps;
	this.precision=precision;
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * function
     */

    private double function(double x){
	return function2use.function(x);
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * derivative of the function
     */

    private double dFunction(double x){
	return function2use.dFunction(x);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * returns the value of the root bracketed between min and max. Starting value of x is determined by 0&lt;=startX&lt;=1
     */

    public double findRoot(double min, double max, double startX, DFunctionOfx function2use, double rightSide){
	int i;
	double deltaX, deltaOld, currentX, aux;
	double f, df, fmin, fmax, result, xdiff;

	this.function2use=function2use;

	fmin=function(min)-rightSide;
	fmax=function(max)-rightSide;
	if(fmin==0) return min;
	if(fmax==0) return max;
	if(fmin*fmax>0){ Output.err.println("Error (in FindRoot/findRoot): Root must be bracketed"); return min; }
	if(fmin>0){
	    aux=min; min=max; max=aux;
	    aux=fmin; fmin=fmax; fmax=aux;
	}
	result=fmax-fmin;
	xdiff=Math.abs(max-min);
	deltaX=deltaOld=xdiff;
	if(startX>1 || startX<0) startX=0.5;
	currentX=min*(1-startX)+max*startX;
	f=function(currentX)-rightSide;
	df=dFunction(currentX);
	for(i=0;i<maxSteps;i++){
	    // Output.err.println("x = "+currentX+" f = "+f+" dx = "+(max-min)+" df = "+(f/df));
	    if(f<0){ min=currentX; fmin=f; }
	    else{ max=currentX; fmax=f; }
	    if(((currentX-max)*df-f)*((currentX-min)*df-f)>0 || Math.abs(2*f)>Math.abs(deltaOld*df)){
		deltaOld=deltaX;
		deltaX=(max-min)/2;
		currentX=min+deltaX;
		if(min==currentX) break;
	    }
	    else{
		deltaOld=deltaX;
		deltaX=f/df;
		aux=currentX;
		currentX-=deltaX;
		if(aux==currentX) break;
	    }
	    if(Math.abs(f)<precision*result) break;
	    // if(Math.abs(df)*xdiff<precision*result){ if(Math.abs(deltaX)<precision*xdiff) break; }
	    // else{ if(Math.abs(df*deltaX)<precision*result) break; }
	    f=function(currentX)-rightSide;
	    df=dFunction(currentX);
	}

	// Output.err.println("Number of steps in FindRoot was "+i);
	// Combine this with "your program" | " awk '/Number of steps in FindRoot/ {a[$(NF)]++} END {for(i in a) print i, a[i]}' |
	// sort -n | awk '{a+=$1*$2; b+=$2; print} END {print "Average is ", a/b}' " to optimize the settings

	if(i==maxSteps) Output.err.println("Warning (in FindRoot/findRoot): Precision "+precision+" has not been reached after "+maxSteps+" steps");
	return currentX;
    }

}
