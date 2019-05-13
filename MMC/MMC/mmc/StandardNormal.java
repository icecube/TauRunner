package mmc;

/**
 * This class provides routines for evaluation of random numbers distributed normally.
 * @author Dmitry Chirkin
 */

public class StandardNormal extends MathModel implements FunctionOfx, FunctionInt{

    private Integral I;
    private static double val1, val2;

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes class with default settings
     */

    public StandardNormal(){
	this(5, 20, 1.e-6);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes the class
     */

    public StandardNormal(int romberg, int maxSteps, double precision){
 	I = new Integral(romberg, maxSteps, precision);
	val1=sndpr(-1);
	val2=sndpr(1);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * evaluates the integrated probability
     */

    public double sndpr(double x){
	if(x<-5) return 0;
	else if(x>5) return 1;
	if(jt) return Math.min(Math.max(J.interpolate(x), 0), 1);
	if(x<=-1) return I.integrateWithSubstitution(1 , x, this, -2);
	else if(x<=1) return val1+I.integrateOpened(-1 , x, this);
	else return val2+I.integrateWithSubstitution(1 , x, this, 2);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * evaluates the standard normal random number
     */

    public double sndrn(double x){
	if(jt) return J.findLimit(x);
	if(x<=val1){
	    I.integrateWithSubstitution(1 , -1, this, -2, -x);
	    return I.getUpperLimit();
	}
	else if(x<=val2){
	    I.integrateOpened(-1 , 1, this, -x+val1);
	    return I.getUpperLimit();
	}
	else{
	    I.integrateWithSubstitution(1 , -1, this, 2, -x+val2);
	    return I.getUpperLimit();
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * evaluates the standard normal random number
     */

    public double sndrn(double rnd, double average, double sigma, double xmin, double xmax, boolean cutoff){
	double x, xl, xh;
	if(xmax<xmin){ x=xmin; xmin=xmax; xmax=x; }
	if(sigma==0) x=average;
	else{
	    if(cutoff) x=rnd;
	    else{
		xl=sndpr((xmin-average)/sigma);
		xh=sndpr((xmax-average)/sigma);
		x=xl+(xh-xl)*rnd;
	    }
	    x=average+sigma*sndrn(x);
	}
	if(x<xmin) x=xmin;
	if(x>xmax) x=xmax;
	return x;
    }

    //----------------------------------------------------------------------------------------------------//

    final private static double norm=1/Math.sqrt(2*Math.PI);

    /**
     * function describes standard normal distribution - interface to Integral
     */

    public double function(double x){
	double aux=norm*Math.exp(-x*x/2);
	return aux;
    }

    //----------------------------------------------------------------------------------------------------//

    public Interpolate J;
    public boolean jt=false;

    /**
     * 1d parametrization - interface to Interpolate
     */

    public double functionInt(double x){
	return sndpr(x);
    }

}
