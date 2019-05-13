package mmc;

/**
 * class contains functions necessary for calculation of photonuclear losses
 */

public class PhotoStochastic extends Photonuclear{

    private Integral[] I;
    private double H[], sum;

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes an integral class separate from that in Propagate
     */

    public PhotoStochastic(Photonuclear cros){
	super(cros);
	int i;
	I = new Integral[m.num];
	for(i=0; i<m.num; i++) I[i] = new Integral(iromb, imaxs, iprec);
	H = new double[m.num];
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * this is what d2N/dv,dx is equal to - interface to Integral
     */

    public double function(double v){
	return cros.cp*photo.photoN(v, cros.component);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * this is what the sum over all components of dN/dx is equal to
     */

    public double dNdx(){
	if(cros.cp<=0) return 0;
	int i;
	sum=0;
	for(i=0; i<m.num; i++){
	    if(jt) sum+=Math.max(Jo[i].interpolate(p.e), 0);
	    else{
		setEnergy(i);
		sum+=I[i].integrateWithLog(vUp, vMax, this);
	    }
	}
	return sum;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * this is what the sum over all components of dN/dx is equal to
     */

    public double dNdx(double rnd){
	if(cros.cp<=0) return 0;
	if(jt) this.rnd=rnd;
	int i;
	sum=0;
	for(i=0; i<m.num; i++){
	    if(jt) H[i]=Math.max(Jo[i].interpolate(p.e), 0);
	    else{
		setEnergy(i);
		H[i]=I[i].integrateWithLog(vUp, vMax, this, rnd);
	    }
	    sum+=H[i];
	}
	return sum;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * this is the value of e=v*E, corresponding to rnd in the call to dNdx
     * for the component chosen with rnd in the argument
     */

    public double e(double rnd){
	int i;
	double rand, rsum;
	rand=rnd*sum;
	rsum=0;
	for(i=0; i<m.num; i++){
	    rsum+=H[i];
	    if(rsum>rand){
		if(jt){
		    setEnergy(i);
		    if(vUp==vMax) return p.e*vUp;
		    return p.e*(vUp*Math.exp(J[i].findLimit(p.e, this.rnd*H[i])*Math.log(vMax/vUp)));
		}
		else{
		    cros.component=i;
		    return p.e*I[i].getUpperLimit();
		}
	    }
	}
	Output.err.println("Error (in PhotoStochastic/e): sum was not initialized correctly");
	return 0;
    }

    //----------------------------------------------------------------------------------------------------//

    public Interpolate J[];
    public boolean jt=false;
    private double rnd;

    /**
     * 2d parametrization - interface to Interpolate
     */

    public double functionInt(double e, double v){
	p.setEnergy(e);
	setEnergy(cros.component);
	if(vUp==vMax) return 0;
	v=vUp*Math.exp(v*Math.log(vMax/vUp));
	return I[cros.component].integrateWithLog(vUp, v, this);
    }

    //----------------------------------------------------------------------------------------------------//

    public Interpolate Jo[];

    /**
     * 1d parametrization - interface to Interpolate
     */

    public double functionInt(double e){
	return J[cros.component].interpolate(e, 1);
    }

}
