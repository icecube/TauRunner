package mmc;

/**
 * class contains functions for evaluation of the spread of the continuous energy losses
 */

public class Energy2LossX extends Energy2Loss{

    private Integral I;

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes an integral class separate from that in Propagate
     */

    public Energy2LossX(Energy2Loss e2loss){
	super(e2loss);
 	I = new Integral(iromb, imaxs, iprec);
    }

    //----------------------------------------------------------------------------------------------------//

    private int name;

    /**
     * total energy^2 losses - interface to Integral
     */

    public double function(double v){
	int i;
	double result;
	switch(name){
	case 0:	return v*v*cros.i.s.function(v);
	case 1: return v*v*cros.b.s.function(v);
	case 2: return v*v*cros.n.s.function(v);
	case 3: return v*v*cros.e.s.function(v);
	default: return 0;
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * contribution of everything to -dE2/dx
     */

    public double dE2dx(){
	if(jt) return Math.max(J.interpolate(p.e), 0);
	int i;
	double sum, aux;
	cros.i.setEnergy(); name=0;
	sum=I.integrateOpened(cros.i.vMin, cros.i.vUp, this);
	for(i=0; i<m.num; i++){
	    cros.b.setEnergy(i); name=1;
	    sum+=I.integrateOpened(0, cros.b.vUp, this);
	}
	for(i=0; i<m.num; i++){
	    cros.n.setEnergy(i); name=2;
	    sum+=I.integrateOpened(cros.n.vMin, cros.n.vUp, this);
	}
	for(i=0; i<m.num; i++){
	    cros.e.setEnergy(i); name=3;
	    sum+=I.integrateOpened(cros.e.vMin, cros.e.vUp, this);
	}
	return p.e*p.e*sum;
    }

    //----------------------------------------------------------------------------------------------------//

    public Interpolate J;
    public boolean jt=false;

    /**
     * 1d parametrization - interface to Interpolate
     */

    public double functionInt(double e){
	p.setEnergy(e);
	return dE2dx();
    }

}
