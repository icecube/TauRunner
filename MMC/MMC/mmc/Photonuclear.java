package mmc;

/**
 * class contains functions necessary for calculation of photonuclear losses
 */

public class Photonuclear extends CrossSections{

    public PhotoContinuous c;
    public PhotoStochastic s;

    public double vMax, vUp, vMin;
    protected Photonuclear photo;

    private Integral I;
    private int i;
    private double v;

    //----------------------------------------------------------------------------------------------------//

    /**
     * creates internal references to p and m, to be called from subclasses
     */

    public Photonuclear(Photonuclear cros){
	p=cros.p;
	m=cros.m;
	this.cros=cros.cros;
	photo=cros;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes subclasses and creates internal references to p and m
     */

    public Photonuclear(CrossSections cros){
	super(cros);
	photo=this;
	c = new PhotoContinuous(this);
	s = new PhotoStochastic(this);
	I = new Integral(iromb, imaxs, iprec);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * call before using the photonuclear functions to set the component of the primary
     */

    public void setEnergy(int i){
	double aux;
	cros.component=i;
	vMin=(Mpi+(Mpi*Mpi)/(2*m.M[i]))/p.e;
	if(p.m<Mpi){
	    aux=p.m/m.M[i];
	    vMax=1-m.M[i]*(1+aux*aux)/(2*p.e);
	}
	else vMax=1;
	vMax=Math.min(vMax, 1-p.m/p.e);
	if(vMax<vMin) vMax=vMin;
	vUp=Math.min(vMax, m.vCut(p.e));
	if(vUp<vMin) vUp=vMin;
	if(photo!=this) photo.vMin=vMin;
    }

    //----------------------------------------------------------------------------------------------------//

    public int bb=1;
    public int form=1;

    /**
     * this is what the photonuclear interaction cross section is equal to
     */

    public double photoN(double v, int i){
	switch(form){
	case 1:
	case 2:
	    {
		final double m1=0.54, m2=1.80;
		double nu, sgn, aux, aum, k, G, t;

		nu=v*p.e;
		nu*=1.e-3;

		switch(bb){
		case 1:
		case 2:
		    if(nu<=200.){
			if(bb==2) sgn=measuredSgN(nu);
			else{
			    if(nu<=17.) sgn=96.1+82./Math.sqrt(nu);
			    else{ aux=Math.log(0.0213*nu); sgn=114.3+1.647*aux*aux; }
			}
		    }
		    else sgn=49.2+11.1*Math.log(nu)+151.8/Math.sqrt(nu);
		    break;
		case 3:  // Bezrukov and Bugaev (BB) parametrization
		    aux=Math.log(0.0213*nu); sgn=114.3+1.647*aux*aux; break;
		case 4:  // ZEUS parametrization
		default:
		    aux=nu*2.e-3*m.M[i]; sgn=63.5*Math.pow(aux, 0.097)+145*Math.pow(aux, -0.5);
		}

		k=1-2/v+2/(v*v);

		if(m.Z[i]==1) G=1;
		else {
		    aux=0.00282*Math.pow(m.A[i], 1./3)*sgn;
		    G=(3/aux)*(0.5+((1+aux)*Math.exp(-aux)-1)/(aux*aux));
		}
		G*=3;

		aux=v*p.m*1.e-3;
		t=aux*aux/(1-v);
		sgn*=1.e-30;

		aum=p.m*1.e-3;
		aum*=aum;
		aux=2*aum/t;
		aux=G*((k+4*aum/m1)*Math.log(1+m1/t)-(k*m1)/(m1+t)-aux)+
		    ((k+2*aum/m2)*Math.log(1+m2/t)-aux)+0.25*aux*((G*(m1-4*t))/(m1+t)+(m2/t)*Math.log(1+t/m2));
		aux*=Alpha/(8*Pi)*m.A[i]*sgn*v;
		if(form==2) if(p.type==1 || p.type==2) aux+=m.A[i]*1.e-30*hardBB(p.e, v);
		return m.No*m.n[i]*p.c*p.c*aux;
	    }
	case 3:
	case 4:
	default:
	    {
		if(jt){
		    setEnergy(i);
		    if(v>=vUp) return Math.max(J[i].interpolate(p.e, Math.log(v/vUp)/Math.log(vMax/vUp)), 0);
		}
		double aux, min, max;
		this.i=i;
		this.v=v;
		min=p.m*v;
		min*=min/(1-v);
		if(p.m<Mpi){
		    aux=p.m*p.m/p.e;
		    min-=(aux*aux)/(2*(1-v));
		}
		max=2*m.M[i]*p.e*(v-vMin);
		//  if(form==4) max=Math.min(max, 5.5e6);  // as requested in Butkevich and Mikheyev
		if(min>max) return 0;
		return m.No*m.n[i]*p.c*p.c*I.integrateWithLog(min, max, this);
	    }
	}
    }

    //----------------------------------------------------------------------------------------------------//

    private boolean initM=true;
    private Interpolate M;

    //----------------------------------------------------------------------------------------------------//

    /*
     * call this to replace photon-nucleon cross section approximation formula with
     * a more exact at low energies parametrization of the experimental data
     */

    private void setMeasured(){
	if(initM){
	    initM=false;
	    final double x[]={0, 0.1, 0.144544, 0.20893, 0.301995, 0.436516, 0.630957, 0.912011, 1.31826, 1.90546, 2.75423, 3.98107, 5.7544, 8.31764, 12.0226, 17.378, 25.1189, 36.3078, 52.4807, 75.8577, 109.648, 158.489, 229.087, 331.131, 478.63, 691.831, 1000, 1445.44, 2089.3, 3019.95, 4365.16, 6309.58, 9120.12, 13182.6, 19054.6, 27542.3, 39810.8, 57544, 83176.4, 120226, 173780, 251188, 363078, 524807, 758576, 1.09648e+06, 1.58489e+06, 2.29086e+06, 3.3113e+06, 4.78628e+06, 6.91828e+06, 9.99996e+06};
	    final double y[]={0, 0.0666667, 0.0963626, 159.74, 508.103, 215.77, 236.403, 201.919, 151.381, 145.407, 132.096, 128.546, 125.046, 121.863, 119.16, 117.022, 115.496, 114.607, 114.368, 114.786, 115.864, 117.606, 120.011, 123.08, 126.815, 131.214, 136.278, 142.007, 148.401, 155.46, 163.185, 171.574, 180.628, 190.348, 200.732, 211.782, 223.497, 235.876, 248.921, 262.631, 277.006, 292.046, 307.751, 324.121, 341.157, 358.857, 377.222, 396.253, 415.948, 436.309, 457.334, 479.025};
	    M = new Interpolate(x, y, 4, false, false);
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * returns measured value of the photon nucleon cross section
     */

    private double measuredSgN(double e){
	setMeasured();
	return M.interpolateArray(e);
    }

    //----------------------------------------------------------------------------------------------------//

    private boolean initH=true;
    private final int hmax=8;
    private Interpolate H[];

    //----------------------------------------------------------------------------------------------------//

    /*
     * initializes hard part of BB (by Bugaev and Shlepin) parametrization calculation
     */

    private void enableHardBB(){
	if(initH){
	    initH=false;
	    final String gcj="";  // gcj would not run with 3d arrays without this
	    final double x[]={3, 4, 5, 6, 7, 8, 9};
	    final double y[][][]={
		{
		    {7.174409e-4, 1.7132e-3, 4.082304e-3, 8.628455e-3, 0.01244159, 0.02204591, 0.03228755},
		    {-0.2436045, -0.5756682, -1.553973, -3.251305, -5.976818, -9.495636, -13.92918},
		    {-0.2942209, -0.68615, -2.004218, -3.999623, -6.855045, -10.05705, -14.37232},
		    {-0.1658391, -0.3825223, -1.207777, -2.33175, -3.88775, -5.636636, -8.418409},
		    {-0.05227727, -0.1196482, -0.4033373, -0.7614046, -1.270677, -1.883845, -2.948277},
		    {-9.328318e-3, -0.02124577, -0.07555636, -0.1402496, -0.2370768, -0.3614146, -0.5819409},
		    {-8.751909e-4, -1.987841e-3, -7.399682e-3, -0.01354059, -0.02325118, -0.03629659, -0.059275},
		    {-3.343145e-5, -7.584046e-5, -2.943396e-4, -5.3155e-4, -9.265136e-4, -1.473118e-3, -2.419946e-3}
		},
		{
		    {-1.269205e-4, -2.843877e-4, -5.761546e-4, -1.195445e-3, -1.317386e-3, -9.689228e-15, -6.4595e-15},
		    {-0.01563032, -0.03589573, -0.07768545, -0.157375, -0.2720009, -0.4186136, -0.8045046},
		    {0.04693954, 0.1162945, 0.3064255, 0.7041273, 1.440518, 2.533355, 3.217832},
		    {0.05338546, 0.130975, 0.3410341, 0.7529364, 1.425927, 2.284968, 2.5487},
		    {0.02240132, 0.05496, 0.144945, 0.3119032, 0.5576727, 0.8360727, 0.8085682},
		    {4.658909e-3, 0.01146659, 0.03090286, 0.06514455, 0.1109868, 0.1589677, 0.1344223},
		    {4.822364e-4, 1.193018e-3, 3.302773e-3, 6.843364e-3, 0.011191, 0.015614, 0.01173827},
		    {1.9837e-5, 4.940182e-5, 1.409573e-4, 2.877909e-4, 4.544877e-4, 6.280818e-4, 4.281932e-4}
		}
	    };
	    H = new Interpolate[hmax];
	    for(int i=0; i<hmax; i++) H[i] = new Interpolate(x, y[p.type-1][i], 4, false, false);
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * returns interpolated value of the hard part of bb cross section
     */

    private double hardBB(double e, double v){
	enableHardBB();
	if(e<1.e5 || v<1.e-7) return 0;
	double aux, sum, lov, loe;
	sum=0; aux=1;
	lov=Math.log(v)/Log10;
	loe=Math.log(e)/Log10-3;
	for(int i=0; i<hmax; i++){
	    if(i>0) aux*=lov;
	    sum+=aux*H[i].interpolateArray(loe);
	}
	return sum/v;
    }

    //----------------------------------------------------------------------------------------------------//

    public int shadow=1;

    /**
     * parametrized photonuclear cross section - interface to Integral
     */

    public double function(double Q2){
	double x, aux, nu, G, F2, R2;
	nu=v*p.e;
	x=Q2/(2*m.M[i]*nu);
	if(m.Z[i]==1) G=1;
	else switch(shadow){
	case 1:
	    {
		if(x<0.0014) G=Math.pow(m.A[i], -0.1);
		else if(x<0.04) G=Math.pow(m.A[i], 0.069*Math.log(x)/Log10+0.097);
		else G=1;
		break;
	    }
	case 2:
	default:
	    {
		if(x>0.3){
		    final double Mb=0.437;
		    final double la=0.5;
		    final double x2=0.278;

		    double mb, Aosc, mu, au, ac;

		    mb=Mb*m.mN[i];
		    au=1/(1-x);
		    ac=1/(1-x2);
		    mu=Mpi/m.M[i];
		    Aosc=(1-la*x)*((au-ac)-mu*(au*au-ac*ac));
		    G=1-mb*Aosc;
		}
		else{
		    final double M1=0.129;
		    final double M2=0.456;
		    final double M3=0.553;

		    double m1, m2, m3, x0, sgn;
		    m1=M1*m.mN[i];
		    m2=M2*m.mN[i];
		    m3=M3*m.mN[i];

		    nu*=1.e-3;
		    sgn=112.2*(0.609*Math.pow(nu, 0.0988)+1.037*Math.pow(nu, -0.5944));

		    aux=0.00282*Math.pow(m.A[i], 1./3)*sgn;
		    G=(3/aux)*(0.5+((1+aux)*Math.exp(-aux)-1)/(aux*aux));
		    G=0.75*G+0.25;

		    x0=Math.pow(G/(1+m2), 1/m1);
		    if(x>=x0) G=Math.pow(x, m1)*(1+m2)*(1-m3*x);
		}
	    }
	}

	switch(form){
	case 3:  // Abramowicz, Levin, Levy, and Maor parametrization
	    {
		double P, W2;

		aux=x*x;
		P=1-1.85*x+2.45*aux-2.35*aux*x+aux*aux;
		G*=(m.Z[i]+(m.A[i]-m.Z[i])*P);

		W2=m.M[i]*m.M[i]-Q2+2*m.M[i]*p.e*v;

		double cp1, cp2, cp3, cr1, cr2, cr3, ap1, ap2, ap3, ar1, ar2, ar3;
		double bp1, bp2, bp3, br1, br2, br3, m2o, m2r, L2, m2p, Q2o;

		switch(bb){
		case 1:  // ALLM91
		    cp1=0.26550;
		    cp2=0.04856;
		    cp3=1.04682;
		    cr1=0.67639;
		    cr2=0.49027;
		    cr3=2.66275;
		    ap1=-0.04503;
		    ap2=-0.36407;
		    ap3=8.17091;
		    ar1=0.60408;
		    ar2=0.17353;
		    ar3=1.61812;
		    bp1=0.49222;
		    bp2=0.52116;
		    bp3=3.55115;
		    br1=1.26066;
		    br2=1.83624;
		    br3=0.81141;
		    m2o=0.30508;
		    m2r=0.20623;
		    L2=0.06527;
		    m2p=10.67564;
		    Q2o=0.27799;
		    break;
		case 2:  // ALLM97
		default:
		    cp1=0.28067;
		    cp2=0.22291;
		    cp3=2.1979;
		    cr1=0.80107;
		    cr2=0.97307;
		    cr3=3.4942;
		    ap1=-0.0808;
		    ap2=-0.44812;
		    ap3=1.1709;
		    ar1=0.58400;
		    ar2=0.37888;
		    ar3=2.6063;
		    bp1=0.60243;
		    bp2=1.3754;
		    bp3=1.8439;
		    br1=0.10711;
		    br2=1.9386;
		    br3=0.49338;
		    m2o=0.31985;
		    m2r=0.15052;
		    L2=0.06527;
		    m2p=49.457;
		    Q2o=0.46017;
		}

		m2o*=1e6;  // GeV -> MeV conversion
		m2r*=1e6;
		L2*=1e6;
		m2p*=1e6;
		Q2o*=1e6;

		bp1*=bp1;  // these values are corrected according to the file f2allm.f from Halina Abramowicz
		bp2*=bp2;
		br1*=br1;
		br2*=br2;
		Q2o+=L2;

		final double R=0;

		double cr, ar, cp, ap, br, bp, t;
		t=Math.log(Math.log((Q2+Q2o)/L2)/Math.log(Q2o/L2));
		if(t<0) t=0;  // gcj 3.x on Linux may create negative values here

		cr=cr1+cr2*Math.pow(t, cr3);
		ar=ar1+ar2*Math.pow(t, ar3);

		cp=cp1+(cp1-cp2)*(1/(1+Math.pow(t, cp3))-1);
		ap=ap1+(ap1-ap2)*(1/(1+Math.pow(t, ap3))-1);

		br=br1+br2*Math.pow(t, br3);
		bp=bp1+bp2*Math.pow(t, bp3);

		double xp, xr, F2p, F2r;

		xp=(Q2+m2p)/(Q2+m2p+W2-m.M[i]*m.M[i]);
		xr=(Q2+m2r)/(Q2+m2r+W2-m.M[i]*m.M[i]);

		F2p=cp*Math.pow(xp, ap)*Math.pow(1-x, bp);
		F2r=cr*Math.pow(xr, ar)*Math.pow(1-x, br);
		F2=(Q2/(Q2+m2o))*(F2p+F2r)*G;
		R2=(2*(1+R));
		break;
	    }
	case 4:  // Butkevich and Mikheyev parametrization
	default:
	    {
		final double a=0.2513e6;
		final double b=0.6186e6;
		final double c=3.0292e6;
		final double d=1.4817e6;
		final double d0=0.0988;
		final double ar=0.4056;
		final double t=1.8152;
		final double As=0.12;
		final double Bu=1.2437;
		final double Bd=0.1853;
		final double R=0.25;

		double F2p, F2n, FSp, FNp, FSn, FNn, n, dl, xUv, xDv;

		n=1.5*(1+Q2/(Q2+c));
		dl=d0*(1+2*Q2/(Q2+d));

		aux=As*Math.pow(x, -dl)*Math.pow(Q2/(Q2+a), 1+dl);
		FSp=aux*Math.pow(1-x, n+4);
		FSn=aux*Math.pow(1-x, n+t);

		aux=Math.pow(x, 1-ar)*Math.pow(1-x, n)*Math.pow(Q2/(Q2+b), ar);
		xUv=Bu*aux;
		xDv=Bd*aux*(1-x);

		FNp=xUv+xDv;
		FNn=xUv/4+xDv*4;

		F2p=FSp+FNp;
		F2n=FSn+FNn;

		F2=G*(m.Z[i]*F2p+(m.A[i]-m.Z[i])*F2n);
		R2=(2*(1+R));
	    }
	}

	aux=Me*Re/Q2;
	aux*=aux*(1-v-m.M[i]*x*v/(2*p.e)+(1-2*p.m*p.m/Q2)*v*v*(1+4*m.M[i]*m.M[i]*x*x/Q2)/R2);
	return (4*Pi*F2/v)*aux;
    }

    //----------------------------------------------------------------------------------------------------//

    public Interpolate J[];
    public boolean jt=false;

    /**
     * 2d parametrization - interface to Interpolate
     */

    public double functionInt(double e, double v){
	p.setEnergy(e);
	setEnergy(cros.component);
	if(vUp==vMax) return 0;
	v=vUp*Math.exp(v*Math.log(vMax/vUp));
	return photoN(v, cros.component);
    }

}
