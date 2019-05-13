package mmc;

/**
 * initializes all cross sections and keeps references to them
 */

public class CrossSections extends PhysicsModel{

    public Decay d;
    public Ionizationloss i;
    public Bremsstrahlung b;
    public Photonuclear n;
    public Epairproduction e;
    public int component=0;

    public double ci=1.;
    public double cb=1.;
    public double cp=1.;
    public double ce=1.;
    public double cd=1.;

    protected Particle p;
    protected Medium m;
    protected CrossSections cros;

    private Integral I;

    //----------------------------------------------------------------------------------------------------//

    /**
     * Necessary to keep class structure, since it is called from subclasses
     */

    public CrossSections(){

    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * creates internal references to p and m, to be called from subclasses
     */

    public CrossSections(CrossSections cros){
	p=cros.p;
	m=cros.m;
	this.cros=cros;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes all cross sections, creates internal references to p and m
     */

    public CrossSections(Particle p, Medium m){
	this.p=p;
	this.m=m;
	d = new Decay(this);
	i = new Ionizationloss(this);
	b = new Bremsstrahlung(this);
	n = new Photonuclear(this);
	e = new Epairproduction(this);
 	I = new Integral(iromb, imaxs, iprec2);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * function for range calculation for given energy - interface to Integral
     */

    public double function(double E){
	final boolean DEBUG=false;
	double result, aux;
	p.setEnergy(E);
	result=0;
	if(DEBUG) Output.err.print(" * "+Output.f(p.e));
	aux=i.c.dEdx();
	result+=aux;
	if(DEBUG) Output.err.print(" \t "+Output.f(aux));
	aux=b.c.dEdx();
	result+=aux;
	if(DEBUG) Output.err.print(" \t "+Output.f(aux));
	aux=n.c.dEdx();
	result+=aux;
	if(DEBUG) Output.err.print(" \t "+Output.f(aux));
	aux=e.c.dEdx();
	result+=aux;
	if(DEBUG) Output.err.println(" \t "+Output.f(aux));
	return -1/result;
    }

    //----------------------------------------------------------------------------------------------------//

    private double ini;

    /**
     * returns the value of the distance integral from ei to ef
     */

    public double getdx(double ei, double ef, double dist){
	if(jt){
	    if(Math.abs(ei-ef)>Math.abs(ei)*halfPrecision){
		ini=J.interpolate(ei);
		double aux=ini-J.interpolate(ef);
		if(Math.abs(aux)>Math.abs(ini)*halfPrecision) return Math.max(aux, 0);
	    }
	    ini=0;
	    return Math.max(Jdf.interpolate((ei+ef)/2)*(ef-ei), 0);
	}
	else return I.integrateWithLog(ei, ef, this, -dist);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * final energy, corresponding to the given value of displacement dist
     */

    public double getef(double ei, double dist){
	if(jt){
	    if(ini!=0){
		double aux=J.findLimit(ini-dist);
		if(Math.abs(aux)>Math.abs(ei)*halfPrecision) return Math.min(Math.max(aux, p.low), ei);
	    }
	    return Math.min(Math.max(ei+dist/Jdf.interpolate(ei+dist/(2*Jdf.interpolate(ei))), p.low), ei);
	}
	else return I.getUpperLimit();
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

    //----------------------------------------------------------------------------------------------------//

    public boolean lpm=false;

}
