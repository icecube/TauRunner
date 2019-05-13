package mmc;

/**
 * class contains functions necessary for calculation of bremsstrahlung losses
 */

public class Bremsstrahlung extends CrossSections{

    public BremsContinuous c;
    public BremsStochastic s;

    public double vMax, vUp, vMin=0;
    protected Bremsstrahlung brems;

    private Integral L;

    //----------------------------------------------------------------------------------------------------//

    /**
     * creates internal references to p and m, to be called from subclasses
     */

    public Bremsstrahlung(Bremsstrahlung cros){
	p=cros.p;
	m=cros.m;
	this.cros=cros.cros;
	brems=cros;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes subclasses and creates internal references to p and m
     */

    public Bremsstrahlung(CrossSections cros){
	super(cros);
	brems=this;
	c = new BremsContinuous(this);
	s = new BremsStochastic(this);
	L = new Integral(iromb, imaxs, iprec);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * call before using the bremsstrahlung functions to set the component of the primary
     */

    public void setEnergy(int i){
	cros.component=i;
	vMax=1-(3./4)*sqrtE*(p.m/p.e)*Math.pow(m.Z[i], 1./3);
	if(vMax<0) vMax=0;
	if(brems.lorenz) vMax=Math.min(vMax, brems.lorenzCut/p.e);
	vMax=Math.min(vMax, 1-p.m/p.e);
	vUp=Math.min(vMax, m.vCut(p.e));
    }

    //----------------------------------------------------------------------------------------------------//

    public int form=1;

    /**
     * this is what the Elastic Bremsstrahlung Cross Section (EBCS) is equal to
     * units are [1/cm] since the multiplication by No*n is done here.
     * Corrections for excitations of the nucleus and deep inelastic excitations of separate nucleons are
     * included (positive term dn/Z), as well as the contribution of the mu-diagrams to the inelastic
     * bremsstrahlung on the electrons (non-zero only for allowed energies of photon after electron recoil).
     */

    protected double Sel(double v, int i){
	double aux, Z3, result, Dn, s1=0;
	Z3=Math.pow(m.Z[i], -1./3);
	switch(form){
	case 1:  // Kelner, Kokoulin, and Petrukhin parametrization
	    {
		int step;
		double d, da, dn, Fa, maxV;
		d=p.m*p.m*v/(2*p.e*(1-v));
		s1=m.B[i]*Z3;
		da=Math.log(1+Me/(d*sqrtE*s1));
		Dn=1.54*Math.pow(m.A[i], 0.27);
		s1=Me*Dn/(p.m*s1);
		dn=Math.log(Dn/(1+d*(Dn*sqrtE-2)/p.m));
		maxV=Me*(p.e-p.m)/(p.e*(p.e-p.p+Me));
		if(v<maxV) Fa=Math.log((p.m/d)/(d*p.m/(Me*Me)+sqrtE))-Math.log(1+Me/(d*sqrtE*m.P[i]*Math.pow(m.Z[i], -2./3)));
		else Fa=0;
		if(m.Z[i]==1) step=0;
		else step=1;
		result=((4./3)*(1-v)+v*v)*(Math.log(p.m/d)-0.5-da-dn+(dn*step+Fa)/m.Z[i]);
	    }
	    break;
	case 2:  // Andreev, Bezrukov, and Bugaev parametrization
	    {
		double aux1, aux2, a1, a2, zeta, qc, qmin, x1, x2, d1, d2, psi1, psi2;
		a1=111.7*Z3/Me; a2=724.2*Z3*Z3/Me;
		qc=1.9*Mmu*Z3;
		aux=2*p.m/qc; aux*=aux;
		zeta=Math.sqrt(1+aux);
		qmin=p.m*p.m*v/(2*p.e*(1-v));
		x1=a1*qmin; x2=a2*qmin;
		if(m.Z[i]==1){ d1=0; d2=0; }
		else{
		    aux1=Math.log(p.m/qc);
		    aux2=(zeta/2)*Math.log((zeta+1)/(zeta-1));
		    d1=aux1+aux2;
		    d2=aux1+((3-zeta*zeta)*aux2+aux)/2;
		}
		aux=p.m*a1; aux1=Math.log(aux*aux/(1+x1*x1));
		aux=p.m*a2; aux2=Math.log(aux*aux/(1+x2*x2));
		psi1=(1+aux1)/2+(1+aux2)/(2*m.Z[i]);
		psi2=(2./3+aux1)/2+(2./3+aux2)/(2*m.Z[i]);
		aux1=x1*Math.atan(1/x1);
		aux2=x2*Math.atan(1/x2);
		psi1-=aux1+aux2/m.Z[i];
		aux=x1*x1; psi2+=2*aux*(1-aux1+(3./4)*Math.log(aux/(1+aux)));
		aux=x2*x2; psi2+=2*aux*(1-aux2+(3./4)*Math.log(aux/(1+aux)))/m.Z[i];
		psi1-=d1; psi2-=d2;
		result=(2-2*v+v*v)*psi1-(2./3)*(1-v)*psi2;
		if(result<0) result=0;
	    }
	    break;
	case 3:  // Petrukhin, Shestakov parametrization
	    {
		double d, Fd;
		d=p.m*p.m*v/(2*p.e*(1-v));
		Fd=189*Z3/Me;
		Fd=p.m*Fd/(1+sqrtE*d*Fd);
		if(m.Z[i]>10) Fd*=(2./3)*Z3;
		result=((4./3)*(1-v)+v*v)*Math.log(Fd);
	    }
	    break;
	case 4:  // complete screening case (Tsai, from PDB)
	default:
	    {
		result=(((4./3)*(1-v)+v*v)*(m.radl(m.Z[i])/m.Z[i])+(1./9)*(1-v)*(m.Z[i]+1))/m.Z[i];
	    }
	}

	aux=2*m.Z[i]*(Me/p.m)*Re;
	aux*=(Alpha/v)*aux*result;
	if(cros.lpm){
	    if(form!=1){
		s1=m.B[i]*Z3;
		Dn=1.54*Math.pow(m.A[i], 0.27);
		s1=Me*Dn/(p.m*s1);
	    }
	    aux*=lpm(v, s1);
	}
	double c2=p.c*p.c;
	return m.No*m.n[i]*c2*c2*aux;
    }

    //----------------------------------------------------------------------------------------------------//

    private boolean init=true;
    private double eLpm;
    protected double Xo;

    /**
     * Landau Pomeranchuk Migdal effect and dielectric suppression evaluation
     */

    protected void setLpm(){
	if(init){
	    boolean lpmSave=cros.lpm;
	    init=false;
	    cros.lpm=false;
	    double sum=0, e=p.e;
	    p.setEnergy(bigEnergy);
	    for(int i=0; i<m.num; i++){
		setEnergy(i);
		sum+=L.integrateOpened(0, vUp, c)+L.integrateWithLog(vUp, vMax, c);
	    }
	    {
		double c2=p.c*p.c;
		Xo=c2*c2/sum;
	    }
	    eLpm=Alpha*p.m;
	    eLpm*=eLpm/(4*Pi*Me*Re*sum);
	    p.setEnergy(e);
	    setEnergy(0);
	    cros.lpm=lpmSave;
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Landau Pomeranchuk Migdal effect and dielectric suppression evaluation
     */

    private double lpm(double v, double s1){
	if(cros.lpm){
	    setLpm();
	    double G, fi, xi, sp, h, s, s2, s3, ps, Gamma;
	    final double fi1=1.54954;
	    final double G1=0.710390;
	    final double G2=0.904912;
	    s1*=s1*sqrt2;
	    sp=Math.sqrt(eLpm*v/(8*p.e*(1-v)));
	    h=Math.log(sp)/Math.log(s1);
	    if(sp<s1) xi=2;
	    else if(sp<1) xi=1+h-0.08*(1-h)*(1-(1-h)*(1-h))/Math.log(s1);
	    else xi=1;
	    s=sp/Math.sqrt(xi);
	    Gamma=Re*Me/(Alpha*p.m*v);
	    Gamma=1+4*Pi*m.No*m.totZ*Re*Gamma*Gamma;
	    s*=Gamma;
	    s2=s*s; s3=s*s2;
	    if(s<fi1) fi=1-Math.exp(-6*s*(1+(3-Pi)*s)+s3/(0.623+0.796*s+0.658*s2));
	    else fi=1-0.012/(s2*s2);
	    if(s<G1){
		ps=1-Math.exp(-4*s-8*s2/(1+3.936*s+4.97*s2-0.05*s3+7.50*s2*s2));
		G=3*ps-2*fi;
	    }
	    else if(s<G2) G=36*s2/(36*s2+1);
	    else G=1-0.022/(s2*s2);
	    return ((xi/3)*((v*v)*G/(Gamma*Gamma)+2*(1+(1-v)*(1-v))*fi/Gamma))/((4./3)*(1-v)+v*v);
	}
	else return 1;
    }

    //----------------------------------------------------------------------------------------------------//

    public boolean lorenz=false;
    public double lorenzCut=1.e6;  // in [MeV]

}
