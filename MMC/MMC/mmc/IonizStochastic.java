package mmc;

/**
 * class contains functions necessary for calculation of ionization losses
 */

public class IonizStochastic extends Ionizationloss{

    private Integral I;

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes an integral class separate from that in Propagate
     */

    public IonizStochastic(Ionizationloss cros){
	super(cros);
 	I = new Integral(iromb, imaxs, iprec);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * this is what d2N/dv,dx is equal to - interface to Integral
     */

    public double function(double v){
	return cros.ci*ioniz.d2Ndvdx(v)*(1+ioniz.inelCorrection(v));
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * this is what dN/dx is equal to
     */

    public double dNdx(){
	if(cros.ci<=0) return 0;
	if(jt) return Math.max(Jo.interpolate(p.e), 0);
	else{
	    setEnergy();
	    return I.integrateWithSubstitution(vUp, vMax, this, 1);
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * this is what dN/dx is equal to
     */

    public double dNdx(double rnd){
	if(cros.ci<=0) return 0;
	if(jt){
	    this.rnd=rnd;
	    return sum=Math.max(Jo.interpolate(p.e), 0);
	}
	else{
	    setEnergy();
	    return I.integrateWithSubstitution(vUp, vMax, this, 1, rnd);
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * this is the value of e=v*p.e, corresponding to rnd in the call to dNdx
     */

    public double e(double rnd){
	int i;
	double rand, rsum;
	rand=m.totZ*rnd;
	rsum=0;
	for(i=0; i<m.num; i++){
	    rsum+=m.n[i]*m.Z[i];
	    if(rsum>rand){
		cros.component=i;
		if(jt){
		    setEnergy();
		    if(vUp==vMax) return p.e*vUp;
		    return p.e*(vUp*Math.exp(J.findLimit(p.e, this.rnd*sum)*Math.log(vMax/vUp)));
		}
		else return p.e*I.getUpperLimit();
	    }
	}
	Output.err.println("Error (in IonizStochastic/e): m.totZ was not initialized correctly");
	return 0;
    }

    //----------------------------------------------------------------------------------------------------//

    public Interpolate J;
    public boolean jt=false;
    private double rnd, sum;

    /**
     * 2d parametrization - interface to Interpolate
     */

    public double functionInt(double e, double v){
	p.setEnergy(e);
	setEnergy();
	if(vUp==vMax) return 0;
	v=vUp*Math.exp(v*Math.log(vMax/vUp));
	return I.integrateWithSubstitution(vUp, v, this, 1);
    }

    //----------------------------------------------------------------------------------------------------//

    public Interpolate Jo;

    /**
     * 1d parametrization - interface to Interpolate
     */

    public double functionInt(double e){
	return J.interpolate(e, 1);
    }

}
