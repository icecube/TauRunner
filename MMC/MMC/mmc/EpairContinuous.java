package mmc;

/**
 * class contains functions necessary for calculation of e+e- pair production losses
 */

public class EpairContinuous extends Epairproduction{

    private Integral I;

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes an integral class separate from that in Propagate
     */

    public EpairContinuous(Epairproduction cros){
	super(cros);
 	I = new Integral(iromb, imaxs, iprec);
    }

    //----------------------------------------------------------------------------------------------------//

    private boolean reverse=false;

    /**
     * pair production energy losses - interface to Integral
     */

    public double function(double v){
	if(reverse) v=1-v;
	return v*epair.ePair(v, cros.component);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * contribution of pair production to -dE/dx
     */

    public double dEdx(){
	if(cros.ce<=0) return 0;
	if(jt) return Math.max(J.interpolate(p.e), 0);
	int i;
	double sum;
	sum=0;
	for(i=0; i<m.num; i++){
	    setEnergy(i);
	    double r1=0.8, rUp=vUp*(1-halfPrecision);
	    boolean rflag=false;
	    if(r1<rUp) if(2*function(r1)<function(rUp)) rflag=true;
	    if(rflag){
		if(r1>vUp) r1=vUp;
		if(r1<vMin) r1=vMin;
		sum+=I.integrateWithLog(vMin, r1, this);
		reverse=true;
		double r2=Math.max(1-vUp, computerPrecision);
		if(r2>1-r1) r2=1-r1;
		sum+=I.integrateOpened(1-vUp, r2, this)+I.integrateWithLog(r2, 1-r1, this);
		reverse=false;
	    }
	    else sum+=I.integrateWithLog(vMin, vUp, this);
	}
	return cros.ce*p.e*sum;
    }

    //----------------------------------------------------------------------------------------------------//

    public Interpolate J;
    public boolean jt=false;

    /**
     * 1d parametrization - interface to Interpolate
     */

    public double functionInt(double e){
	p.setEnergy(e);
	return dEdx();
    }

}
