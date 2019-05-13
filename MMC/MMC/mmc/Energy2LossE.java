package mmc;

/**
 * class contains functions for evaluation of the spread of the continuous energy losses
 */

public class Energy2LossE extends Energy2Loss{

    private Integral I;

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes an integral class separate from that in Propagate
     */

    public Energy2LossE(Energy2Loss e2loss){
	super(e2loss);
 	I = new Integral(iromb, imaxs, iprec2);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * total energy^2 losses - interface to Integral
     */

    public double function(double E){
	double aux;
	aux=cros.function(E);
	return aux*e2loss.e2lx.dE2dx();
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * contribution of everything to -dE2/de
     */

    public double dE2de(double ei, double ef){
	if(jt){
	    if(Math.abs(ei-ef)>Math.abs(ei)*halfPrecision){
		double aux=J.interpolate(ei);
		double aux2=aux-J.interpolate(ef);
		if(Math.abs(aux2)>Math.abs(aux)*halfPrecision) return Math.max(aux2, 0);
	    }
	    else return Math.max(Jdf.interpolate((ei+ef)/2)*(ef-ei), 0);
	}
	return I.integrateWithLog(ei, ef, this);
    }

    //----------------------------------------------------------------------------------------------------//

    boolean df=false;
    public Interpolate J;
    public Interpolate Jdf;
    public boolean jt=false;

    /**
     * 1d parametrization - interface to Interpolate
     */

    public double functionInt(double e){
	if(df) return function(e);
	else return I.integrateWithLog(e, p.low, this);
    }

}
