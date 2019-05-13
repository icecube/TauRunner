package mmc;

/**
 * particle definition (muon)
 */

public class Particle extends PhysicsModel{

    public double r=0, x=0, y=0, z=0, t=0;             // coordinates [cm] and age [sec]
    public double theta=0, phi=0;                      // zenith and azimuth of the momentum in [deg]
    public double costh=1, sinth=0, cosph=1, sinph=0;  // cos and sin of theta and phi
    public double p=0, p2=0, e=0;                      // momentum, momentum square and energy in [MeV]
    public double m=Mmu, l=Lmu, c=1;                   // mass [MeV], lifetime [sec] and charge
    public String name="mu";                           // name of the particle
    public double low=m;                               // energy below which the particle is lost [MeV]
    public int type=1;                                 // particle type: 1 for muon, 2 for tau, 3 for electron
    public int igen=0, gens=1;                         // parent particle id, particle id

    protected Scattering s;
    protected Propagate pr;

    private Integral I;

    public double xi=0, yi=0, zi=0, ti=0, Ei=0;        // coordinates at entry point [m,m,m,sec,GeV]
    public double xf=0, yf=0, zf=0, tf=0, Ef=0;        // coordinates at exit point
    public double xc=0, yc=0, zc=0, tc=0, Ec=0;        // coordinates at point of closest approach
    public double Elost=0;                             // energy lost in the detector volume [GeV]

    //----------------------------------------------------------------------------------------------------//

    /**
     * initialize particle
     */

    public Particle(Propagate pr, String name){
	if(name.startsWith("tau")){
	    this.name="tau";
	    type=2;
	    m=Mtau;
	    l=Ltau;
	}
	else if(name.startsWith("e")){
	    this.name="e";
	    type=3;
	    m=Me;
	    l=-1;
	}
	else if(name.startsWith("monopole")){
	    type=4;
	    try{
		m=Double.parseDouble(name.substring(9))*1.e3;
	    }catch(Exception error){
		m=0;
	    }
	    if(m<=0) m=Mmon;
	    this.name="monopole-"+Output.f(m*1.e-3);
	    c=Cmon;
	    l=-1;
	}
	else if(name.startsWith("stau")){
	    type=5;
	    try{
		m=Double.parseDouble(name.substring(5))*1.e3;
	    }catch(Exception error){
		m=0;
	    }
	    if(m<=0) m=Mstau;
	    this.name="stau-"+Output.f(m*1.e-3);
	    l=Lstau;
	}
	else{
	    this.name="mu";
	    type=1;
	    m=Mmu;
	    l=Lmu;
	}
	low=m; e=m;
	this.pr=pr;
	if(elow>low) low=elow;
	if(ebig<100*low) ebig=Math.max(100*low, bigEnergy);
	s = new Scattering(this);
 	I = new Integral(iromb, imaxs, iprec2);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * store particle information
     */

    public Particle(int igen, int gens, String name, double x, double y, double z, double theta, double phi, double e, double t, double r, Particle p){
	this(name, x, y, z, theta, phi, e, t);
	this.r=r;
	this.igen=igen;
	this.gens=gens;
	if(p!=null){
	    xi=p.xi; yi=p.yi; zi=p.zi; ti=p.ti; Ei=p.Ei;
	    xf=p.xf; yf=p.yf; zf=p.zf; tf=p.tf; Ef=p.Ef;
	    xc=p.xc; yc=p.yc; zc=p.zc; tc=p.tc; Ec=p.Ec;
	    Elost=p.Elost;
	}
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * store particle information
     */

    public Particle(int igen, int gens, String name, double x, double y, double z, double theta, double phi, double e, double t, double r){
	this(name, x, y, z, theta, phi, e, t);
	this.r=r;
	this.igen=igen;
	this.gens=gens;
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * store particle information
     */

    public Particle(String aname, double x, double y, double z, double theta, double phi, double e, double t){
	String name=aname.length()==0?"?":aname.charAt(0)=='a'?aname.substring(1):aname;
	if(name.equals("tau") || name.equals("tau-") || name.equals("tau+")){
	    if(name.equals("tau+")) type=-33;
	    else type=-34;
	    m=Mtau;
	    l=Ltau;
	}
	else if(name.equals("mu") || name.equals("mu-") || name.equals("mu+")){
	    if(name.equals("mu+")) type=-5;
	    else type=-6;
	    m=Mmu;
	    l=Lmu;
	}
	else if(name.startsWith("stau") ){
	    if(name.indexOf("stau+")!=-1) type=-9131; else type=-9132;
	    try{
		m=Double.parseDouble(name.substring(5))*1.e3;
	    }catch(Exception error){
		m=0;
	    }
	    if(m<=0) m=Mstau;
	    l=Lstau;
	}
	else if(name.equals("e") || name.equals("e-") || name.equals("e+")){
	    if(name.equals("e+")) type=-2;
	    else type=-3;
	    m=Me;
	    l=-1;
	}
	else if(name.indexOf("nu_")!=-1){
	    if(name.equals("nu_e")) type=-201;
	    else if(name.equals("~nu_e")) type=-204;
	    else if(name.equals("nu_mu")) type=-202;
	    else if(name.equals("~nu_mu")) type=-205;
	    else if(name.equals("nu_tau")) type=-203;
	    else if(name.equals("~nu_tau")) type=-206;
	    m=0;
	    l=-1;
	}
	else{
	    if(name.equals("delta")) type=-1002;
	    else if(name.equals("brems")) type=-1001;
	    else if(name.equals("munu")) type=-1004;
	    else if(name.equals("epair")) type=-1003;
	    else if(name.equals("hadr")) type=-1006;
	    else if(name.equals("conti")) type=-1111;
	    else type=0;
	    m=0;
	    l=0;
	}
	low=m;
	setEnergy(e);
	location(aname, t, x, y, z, theta, phi);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initialize the location and direction of the particle, time in sec, x, y, z in cm, theta and phi in deg
     */

    public void location(String name, double time, double x, double y, double z, double theta, double phi){
	this.name=name;
	r=0;
	t=time;
	this.x=x;
	this.y=y;
	this.z=z;
	this.theta=theta;
	this.phi=phi;
	theta*=(Pi/180);
	phi*=(Pi/180);
	costh=Math.cos(theta);
	sinth=Math.sin(theta);
	cosph=Math.cos(phi);
	sinph=Math.sin(phi);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * advances the particle by the given distance
     */

    public void advance(double dr, double ei, double ef){
	r+=dr;
	if(pr.exactTime) t+=getdt(ei, ef)/pr.rho;
	else t+=dr/C;
	if(pr.molieScat){
	    double tho, max, rnd1, rnd2, sx, sy, sz, tx, ty, tz;
	    tho=s.gettho(dr, ei, ef);
	    max=1/sqrt2;

	    rnd1=pr.S.sndrn(Propagate.rand.nextDouble(), 0, tho, -max, max, false);
	    rnd2=pr.S.sndrn(Propagate.rand.nextDouble(), 0, tho, -max, max, false);

	    sx=(rnd1/sqrt3+rnd2)/2;
	    tx=rnd2;

	    rnd1=pr.S.sndrn(Propagate.rand.nextDouble(), 0, tho, -max, max, false);
	    rnd2=pr.S.sndrn(Propagate.rand.nextDouble(), 0, tho, -max, max, false);

	    sy=(rnd1/sqrt3+rnd2)/2;
	    ty=rnd2;

	    sz=Math.sqrt(Math.max(1-(sx*sx+sy*sy), 0));
	    tz=Math.sqrt(Math.max(1-(tx*tx+ty*ty), 0));

	    double ax, ay, az;

	    ax=sinth*cosph*sz+costh*cosph*sx-sinph*sy;
	    ay=sinth*sinph*sz+costh*sinph*sx+cosph*sy;
	    az=costh*sz-sinth*sx;

	    x+=ax*dr;
	    y+=ay*dr;
	    z+=az*dr;

	    ax=sinth*cosph*tz+costh*cosph*tx-sinph*ty;
	    ay=sinth*sinph*tz+costh*sinph*tx+cosph*ty;
	    az=costh*tz-sinth*tx;

	    costh=az;
	    sinth=Math.sqrt(Math.max(1-costh*costh, 0));
	    if(sinth!=0){
		sinph=ay/sinth;
		cosph=ax/sinth;
	    }

	    theta=Math.acos(costh>1?1:costh<-1?-1:costh)*180/Pi;
	    phi=Math.acos(cosph>1?1:cosph<-1?-1:cosph)*180/Pi;
	    if(sinph<0) phi=360-phi; if(phi>=360) phi-=360;

	}
	else{
	    x+=sinth*cosph*dr;
	    y+=sinth*sinph*dr;
	    z+=costh*dr;
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * sets the energy of the particle
     */

    public void setEnergy(double e){
	this.e=e;
	p2=e*e-m*m;
	p=Math.sqrt(Math.max(p2, 0));
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * function for time delta calculation - interface to Integral
     */

    public double function(double E){
	final boolean DEBUG=false;
	double aux;
	aux=pr.s.function(E);
	aux*=e/(p*C);
	if(DEBUG) Output.err.println(" # "+Output.f(e)+" \t "+Output.f(aux));
	return aux;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * time delta, corresponding to the given propagation distance
     */

    public double getdt(double ei, double ef){
	if(jt){
	    if(Math.abs(ei-ef)>Math.abs(ei)*halfPrecision){
		double aux=J.interpolate(ei);
		double aux2=aux-J.interpolate(ef);
		if(Math.abs(aux2)>Math.abs(aux)*halfPrecision) return aux2;
	    }
	    return Jdf.interpolate((ei+ef)/2)*(ef-ei);
	}
	else return I.integrateWithLog(ei, ef, this);
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
	else return I.integrateWithLog(e, low, this);
    }

}
