package mmc;

/**
 * class contains functions for evaluation of the Moliere angle distribution spread
 */

public class Scattering extends PhysicsModel{

    private Particle p;
    private Propagate pr;
    private Integral I;

    //----------------------------------------------------------------------------------------------------//

    /**
     * initialize the class
     */

    public Scattering(Particle p){
	this.p=p;
	pr=p.pr;
 	I = new Integral(iromb, imaxs, iprec2);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * function for angle distribution spread calculation - interface to Integral
     */

    public double function(double E){
	final boolean DEBUG=false;
	double aux, aux2;
	aux=pr.s.function(E);
	aux2=Ry*p.e/(p.p*p.p);
	aux*=aux2*aux2;
	if(DEBUG) Output.err.println(" $ "+Output.f(p.e)+" \t "+Output.f(aux));
	return aux;
    }

    //----------------------------------------------------------------------------------------------------//

    final private static double cutoff=1;

    /**
     * spread of the angle distribution, corresponding to the given propagation distance
     */

    public double gettho(double dr, double ei, double ef){
	double aux;
	if(jt){
	    if(Math.abs(ei-ef)>Math.abs(ei)*halfPrecision){
		aux=J.interpolate(ei);
		double aux2=aux-J.interpolate(ef);
		if(Math.abs(aux2)>Math.abs(aux)*halfPrecision) aux=aux2;
		else aux=Jdf.interpolate((ei+ef)/2)*(ef-ei);
	    }
	    else aux=Jdf.interpolate((ei+ef)/2)*(ef-ei);
	}
	else aux=I.integrateWithLog(ei, ef, this);
	double Xo=pr.m.Xo;
	// pr.s.b.setLpm(); Xo=pr.s.b.Xo; // correct only for electrons with CSC
	aux=Math.sqrt(Math.max(aux, 0)/Xo)*p.c*Math.max(1+0.038*Math.log(dr/Xo), 0);
	return Math.min(aux, cutoff);
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
	else return I.integrateWithLog(e, ebig, this);
    }

}
