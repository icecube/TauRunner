package mmc;

/**
 * class contains functions necessary for the calculation of decay
 */

public class Decay extends CrossSections implements DFunctionOfx{

    private FindRoot f;

    /**
     * creates internal references to p and m
     */

    public Decay(CrossSections cros){
	super(cros);
	f = new FindRoot(imaxs, iprec);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * this cross section describes decay
     */

    public double decay(){
	if(cros.cd<=0 || p.l<0) return 0;
	return cros.cd/Math.max((p.p/p.m)*p.l*C, xres);
    }

    //----------------------------------------------------------------------------------------------------//

    public static boolean flag=false;
    public String out;

    /**
     * energy of the electron that results from the muon decay
     */

    public double e(double ernd, double arnd, double srnd, Output o){
	if(p.l<0) return 0;
	double emax, x0, f0, el, lm, pl;
	String out1="nu", out2="nu";
	if(p.type==2){
	    final double brmu=0.1737;
	    final double brel=0.1783+brmu;
	    final double brpi=0.1109+brel;
	    final double br2p=0.2540+brpi;
	    final double br3p=0.1826+br2p;
	    if(srnd<brmu){
		lm=Mmu;
		if(p.name.endsWith("+")){ out="mu+"; if(flag){ out1="~nu_tau"; out2="nu_mu"; } }
		else if(p.name.endsWith("-")){ out="mu-"; if(flag){ out1="nu_tau"; out2="~nu_mu"; } }
		else{ out="mu"; if(flag){ out1="nu_tau"; out2="~nu_mu"; } }
	    }
	    else if(srnd<brel){
		lm=Me;
		if(p.name.endsWith("+")){ out=Output.AMASIM?"epair":"e+"; if(flag){ out1="~nu_tau"; out2="nu_e"; } }
		else if(p.name.endsWith("-")){ out=Output.AMASIM?"delta":"e-"; if(flag){ out1="nu_tau"; out2="~nu_e"; } }
		else{ out=Output.AMASIM?"delta":"e"; if(flag){ out1="nu_tau"; out2="~nu_e"; } }
	    }
	    else{
		if(srnd<brpi) lm=Mpi;
		else if(srnd<br2p) lm=Mrh;
		else if(srnd<br3p) lm=Ma1;
		else lm=Mrs;
		el=(Mtau*Mtau+lm*lm)/(2*Mtau);
		out=Output.AMASIM?"munu":"hadr";
		el=el*(p.e/p.m)+Math.sqrt(el*el-lm*lm)*(p.p/p.m)*(2*arnd-1);
		if(flag){
		    if(p.name.endsWith("+")){ out1="~nu_tau"; }
		    else if(p.name.endsWith("-")){ out1="nu_tau"; }
		    else{ out1="nu_tau"; }
		    o.output(1, out1, p.e-el, 0);
		}
		return el;
	    }
	}
	else{
	    lm=Me;
	    if(p.name.endsWith("+")){ out=Output.AMASIM?"epair":"e+"; if(flag){ out1="~nu_mu"; out2="nu_e"; } }
	    else if(p.name.endsWith("-")){ out=Output.AMASIM?"delta":"e-"; if(flag){ out1="nu_mu"; out2="~nu_e"; } }
	    else{ out=Output.AMASIM?"delta":"e"; if(flag){ out1="nu_mu"; out2="~nu_e"; } }
	}
	emax=(p.m*p.m+lm*lm)/(2*p.m);
	x0=lm/emax;
	f0=x0*x0;
	f0=f0*x0-f0*f0/2;
	el=Math.max(f.findRoot(x0, 1, 0.5, this, f0+(0.5-f0)*ernd)*emax, lm);
	pl=Math.sqrt(el*el-lm*lm);
	if(flag){
	    int sign;
	    double cp, sp, ct, st, en, En1, En2;
	    cp=pl/(p.m-el);
	    sp=Math.sqrt(1-cp*cp);
	    if(arnd<0.5){
		sign=1;
		arnd=2*arnd;
	    }
	    else{
		sign=-1;
		arnd=2*arnd-1;
	    }
	    ct=2*arnd-1;
	    st=Math.sqrt(1-ct*ct);
	    en=(p.m-el)/2;
	    En1=-cp*ct-sign*sp*st;
	    En2=-cp*ct+sign*sp*st;
	    En1=en*((p.e/p.m)+(p.p/p.m)*En1);
	    En2=en*((p.e/p.m)+(p.p/p.m)*En2);
	    o.output(1, out1, En1, 0);
	    o.output(1, out2, En2, 0);
	}
	return el*(p.e/p.m)+pl*(p.p/p.m)*(2*arnd-1);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * function for electron energy calculation - interface to FindRoot
     */

    public double function(double x){
	double x2;
	x2=x*x;
	return x*x2-x2*x2/2;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * function for electron energy calculation - interface to FindRoot
     */

    public double dFunction(double x){
	return (3-2*x)*x*x;
    }

}
