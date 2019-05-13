package tfa;
import mmc.*;
import java.io.Reader;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.util.Vector;
import java.util.StringTokenizer;

/**
 * Implements muon/tau propagation through multiple media to/through the detector.
 * Current possibilities for detector geometry are cylinder, cuboid, and sphere.
 * The origin is at the center of the detector.
 */

public class Amanda{

    private final static boolean lfix=true;
    private final static double LENGTH=800;
    private final static double RADIUS=400;
    private final static double WIDTH=800;
    private final static double HEIGHT=800;
    private final static double BIG=6.40e7;

    private boolean SURF=false;
    private boolean FACE=false;
    private boolean USER=false;
    private boolean USFI=false;
    private boolean SDEC=false;
    private boolean RECC=false;
    private double zset=0;
    private String type;

    private int mediamax=100;
    private String mediadef=null;

    private double vcut[] = { -1.0, -1.0 };
    private double ecut[] = { -1.0, -1.0 };
    private String med="Ice";
    private String muta="mu";
    private String usna="mmc_en";
    private double elow=PhysicsModel.elow;
    private double ebig=PhysicsModel.ebig;

    private boolean conti[] = { false, false };
    private boolean timef=false;
    private boolean lpmef=false;
    private boolean scatt=false;
    private boolean frho=false;
    private boolean rfix=false;
    private boolean amasim=false;

    private int bspar=1;
    private int pncrs=1;
    private int pncbb=1;
    private int pncsh=1;

    private double crsci=1;
    private double crscb=1;
    private double crscp=1;
    private double crsce=1;
    private double crscd=1;
    private double drho=1;

    private int romb=5;
    private long seed=0;
    private boolean SEED=false;
    private boolean raw=false;
    private String tdir="";
    public String dtct="";

    public int gdet=0;
    public double length=LENGTH, radius=RADIUS, width=WIDTH, height=HEIGHT;
    public double rho[], srho[];
    private Propagate p1[], p2[], p3[];
    private double surf=0;
    private double emax=0;
    private StringBuffer hist;
    private Vector userbf;
    private Vector I3hist;
    private Particle pI;

    private int igen, gens=0, imax=0;
    private String prnt;
    public long events=0, tracks=0, vertices=0, missed=0;

    private String param="Amanda", options;
    private int fnu;
    private int fnm[];
    private double fnx[];
    public int medt[];
    public double sphz[];
    public double sphr[];
    public double boxx[];
    public double boxy[];
    public double boxz[];
    public double boxl[];
    public double boxw[];
    public double boxh[];
    public double cylz[];
    public double cylr[];
    public double cyll[];
    public int medianum=1;
    public boolean DEBUG=false;

    //----------------------------------------------------------------------------------------------------//

    /**
     * Front-end for AMANDA, reads from and outputs to standard input/output in F2000 format.
     */

    public static void main(String[] args){
	Amanda A = new Amanda();
	A.mmcf2k(args);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initializes propagation for external applications, e.g.&nbsp;mmc-icetray (through jni).
     */

    public void setup(String args){
	dtct=".icecube";
	setup(Output.splitString(args));
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initializes propagation for external applications, e.g.&nbsp;mmc-icetray (through jni).
     */

    public void setup(String args[]){
	Output.I3flag=true;
	setp(args);
	I3hist = new Vector(Output.HISTSIZE);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Propagates particles for the external applications. Returns an array of secondaries.
     * If "-user" option is used, fills user-block variables.
     */

    public Particle[] propagate(Particle p){
	if(!Output.I3flag) return null;
	Particle[] I3p;

	pI=p;
	type=p.name;
	if(muta.equals(type) || (muta+"-").equals(type) || (muta+"+").equals(type)){
	    p.r=prop(p.igen, p.gens, p.x*1.e-2, p.y*1.e-2, p.z*1.e-2, 180-p.theta, p.phi<180?p.phi+180:p.phi-180, p.r*1.e-2, p.e*1.e-3, p.t*1.e9);
	    I3p = new Particle[I3hist.size()];
	    for(int i=0; i<I3hist.size(); i++) I3p[i]=(Particle)I3hist.elementAt(i);
	}
	else{
	    I3p = new Particle[0];
	}
	return I3p;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * This is the command-line option parser. Call with "-help" to list all options.
     */

    public boolean setp(String[] args){
	int i, j;

	int bnum=0;
	String pbad="";
	boolean pflag;

	rho = new double[mediamax];
	srho = new double[mediamax];
	String med[] = new String[mediamax];
	double vcut[][] = { new double[mediamax], new double[mediamax] };
	double ecut[][] = { new double[mediamax], new double[mediamax] };
	boolean conti[][] = { new boolean[mediamax], new boolean[mediamax] };

	rho[0]=this.drho;
	med[0]=this.med;
	vcut[0][0]=this.vcut[0];
	vcut[1][0]=this.vcut[1];
	ecut[0][0]=this.ecut[0];
	ecut[1][0]=this.ecut[1];
	conti[0][0]=this.conti[0];
	conti[1][0]=this.conti[1];

	String mediaopts="";
	int mediaseg=1;
	medt = new int[mediamax];
	sphz = new double[mediamax];
	sphr = new double[mediamax];
	boxx = new double[mediamax];
	boxy = new double[mediamax];
	boxz = new double[mediamax];
	boxl = new double[mediamax];
	boxw = new double[mediamax];
	boxh = new double[mediamax];
	cylz = new double[mediamax];
	cylr = new double[mediamax];	
	cyll = new double[mediamax];	
	boolean medn[] = new boolean[mediamax];

	medianum=1;
	medt[0]=0;
	medn[0]=true;

	for(int n=0; n<args.length; n++){
	    pflag=true;
	    if(args[n].equals("-help") || args[n].equals("-h") || args[n].equals("--help")){
		Output.out.println("\n"+
"This program propagates muons in "+med[0]+" to/through the detector\n"+
"Available options are: -length=[LENGTH of the detector volume in meters]\n"+
"                       -radius=[RADIUS of the detector volume in meters]\n"+
"                       -width=[WIDTH of the detector volume in meters]\n"+
"                       -height=[HEIGHT of the detector volume in meters]\n"+
"                       -gdet=[0-2] detector is a cylinder/box/sphere\n"+
"                       -vcut=[value of vcut used for the 1st region]\n"+
"                       -ecut=[ecut in  MeV  used for the 2nd region]\n"+
"                       -medi=[medium name]\n"+
"                       -mediadef=[file with media definitions]\n"+
"                       -tau  propagate taus instead of muons\n"+
"                       -e    propagate electrons instead of muons\n"+
"                       -monopole[=mass in GeV] propagate monopoles\n"+
"                       -stau[=mass in GeV]     propagate staus\n"+
"                       -sdec enable stopped muon decay treatment\n"+
"                       -recc enable printout of continuous energy losses\n"+
"                       -user     enable the mmc_en user line\n"+
"                       -user=[z] same, but record energy at z, not CPD\n"+
"                       -rdmc     enforce compliance with rdmc\n"+
"                       -amasim   turn on workarounds for amasim\n"+
"                       -time precise time of flight calculation\n"+
"                       -cont enable continuous loss randomization\n"+
"                       -scat enable Moliere scattering\n"+
"                       -lpm  enable lpm treatment\n"+
"                       -bs=[1-4]  bremsstrahlung: kkp, abb, ps, csc\n"+
"                       -ph=[1-4]  photonuclear: bb, bb+bs, allm, bm\n"+
"                       -bb=[bb/bs:1-4 3|4:bb|zeus, allm:1-2(91/7), bm:1]\n"+
"                       -sh=[1-2]  nuclear structure function: dutt/butk\n"+
"                       -c[i/b/p/e/d]=[cross section modifier, 0:disable]\n"+
"                       -rho=[multiplicative factor] for medium density\n"+
"                       -frho  enable smart density factor handling\n"+
"                       -rw=[reweight cross sections by 1-x^+rw or x^-rw]\n"+
"                       -elow=[muon energy in GeV below which it is lost]\n"+
"                       -ebig=[upper bound in GeV of the paramet. tables]\n"+
"                       -surf=[h in meters]  propagate to the plane z=[h]\n"+
"                       -face  only if detector is on opposite side of it\n"+
"                       -romb=[number of interpolation points]\n"+
"                       -seed=[integer] sets random number generator seed\n"+
"                       -raw  save tables in raw format\n"+
"                       -tdir=[dir] specify directory for paramet. tables\n");
		return false;
	    }
	    else if(args[n].equals("-tau")){
		muta="tau";
		usna="mmc_et";
	    }
	    else if(args[n].equals("-e")){
		muta="e";
		usna="mmc_el";
	    }
	    else if(args[n].equals("-sdec")){
		SDEC=true;
	    }
	    else if(args[n].equals("-recc")){
		RECC=true;
	    }
	    else if(args[n].equals("-user")){
		USER=true;
	    }
	    else if(args[n].equals("-rdmc")){
		rfix=true;
	    }
	    else if(args[n].equals("-amasim")){
		amasim=true;
	    }
	    else if(args[n].equals("-time")){
		timef=true;
	    }
	    else if(args[n].equals("-cont")){
		conti[0][0]=true;
	    }
	    else if(args[n].equals("-scat")){
		scatt=true;
	    }
	    else if(args[n].equals("-lpm")){
		lpmef=true;
	    }
	    else if(args[n].equals("-frho")){
		frho=true;
	    }
	    else if(args[n].equals("-face")){
		FACE=true;
	    }
	    else if(args[n].equals("-raw")){
		raw=true;
	    }
	    else if(args[n].startsWith("-monopole")){
		muta=args[n].substring(1);
		usna="mmc_mn";
	    }
	    else if(args[n].startsWith("-stau")){
		muta=args[n].substring(1);
		usna="mmc_st";
	    }
	    else if(args[n].startsWith("-length=")){
		try{
		    length=Double.parseDouble(args[n].substring(8));
		}catch(NumberFormatException error){
		    length=0;
		}
	    }
	    else if(args[n].startsWith("-radius=")){
		try{
		    radius=Double.parseDouble(args[n].substring(8));
		}catch(NumberFormatException error){
		    radius=0;
		}
	    }
	    else if(args[n].startsWith("-width=")){
		try{
		    width=Double.parseDouble(args[n].substring(7));
		}catch(NumberFormatException error){
		    width=0;
		}
	    }
	    else if(args[n].startsWith("-height=")){
		try{
		    height=Double.parseDouble(args[n].substring(8));
		}catch(NumberFormatException error){
		    height=0;
		}
	    }
	    else if(args[n].startsWith("-surf=")){
		try{
		    surf=Double.parseDouble(args[n].substring(6));
		}catch(NumberFormatException error){
		    surf=0;
		}
		SURF=true;
	    }
	    else if(args[n].startsWith("-vcut=")){
		try{
		    double aux=Double.parseDouble(args[n].substring(6));
		    if(aux<0) ecut[0][0]=-aux; else vcut[0][0]=aux;
		}catch(NumberFormatException error){
		    vcut[0][0]=this.vcut[0];
		}
	    }
	    else if(args[n].startsWith("-ecut=")){
		try{
		    double aux=Double.parseDouble(args[n].substring(6));
		    if(aux<0) vcut[1][0]=-aux; else ecut[1][0]=aux;
		}catch(NumberFormatException error){
		    ecut[1][0]=this.ecut[1];
		}
	    }
	    else if(args[n].startsWith("-user=")){
		try{
		    USER=true; USFI=true;
		    zset=Double.parseDouble(args[n].substring(6));
		}catch(NumberFormatException error){
		    zset=0;
		}
	    }
	    else if(args[n].startsWith("-gdet=")){
		try{
		    gdet=(int)Double.parseDouble(args[n].substring(6));
		}catch(NumberFormatException error){
		    gdet=0;
		}
	    }
	    else if(args[n].startsWith("-bs=")){
		try{
		    bspar=(int)Double.parseDouble(args[n].substring(4));
		}catch(NumberFormatException error){
		    bspar=1;
		}
	    }
	    else if(args[n].startsWith("-ph=")){
		try{
		    pncrs=(int)Double.parseDouble(args[n].substring(4));
		}catch(NumberFormatException error){
		    pncrs=1;
		}
	    }
	    else if(args[n].startsWith("-bb=")){
		try{
		    pncbb=(int)Double.parseDouble(args[n].substring(4));
		}catch(NumberFormatException error){
		    pncbb=1;
		}
	    }
	    else if(args[n].startsWith("-sh=")){
		try{
		    pncsh=(int)Double.parseDouble(args[n].substring(4));
		}catch(NumberFormatException error){
		    pncsh=1;
		}
	    }
	    else if(args[n].startsWith("-ci=")){
		try{
		    crsci=Double.parseDouble(args[n].substring(4));
		}catch(NumberFormatException error){
		    crsci=1;
		}
	    }
	    else if(args[n].startsWith("-cb=")){
		try{
		    crscb=Double.parseDouble(args[n].substring(4));
		}catch(NumberFormatException error){
		    crscb=1;
		}
	    }
	    else if(args[n].startsWith("-cp=")){
		try{
		    crscp=Double.parseDouble(args[n].substring(4));
		}catch(NumberFormatException error){
		    crscp=1;
		}
	    }
	    else if(args[n].startsWith("-ce=")){
		try{
		    crsce=Double.parseDouble(args[n].substring(4));
		}catch(NumberFormatException error){
		    crsce=1;
		}
	    }
	    else if(args[n].startsWith("-cd=")){
		try{
		    crscd=Double.parseDouble(args[n].substring(4));
		}catch(NumberFormatException error){
		    crscd=1;
		}
	    }
	    else if(args[n].startsWith("-rho=")){
		try{
		    rho[0]=Double.parseDouble(args[n].substring(5));
		}catch(NumberFormatException error){
		    rho[0]=this.drho;
		}
	    }
	    else if(args[n].startsWith("-rw=")){
		try{
		    rw=Double.parseDouble(args[n].substring(4));
		}catch(NumberFormatException error){
		    rw=0;
		}
	    }
	    else if(args[n].startsWith("-elow=")){
		try{
		    elow=Double.parseDouble(args[n].substring(6))*1.e3;
		}catch(NumberFormatException error){
		    elow=PhysicsModel.elow;
		}
	    }
	    else if(args[n].startsWith("-ebig=")){
		try{
		    ebig=Double.parseDouble(args[n].substring(6))*1.e3;
		}catch(NumberFormatException error){
		    ebig=PhysicsModel.ebig;
		}
	    }
	    else if(args[n].startsWith("-romb=")){
		try{
		    romb=(int)Double.parseDouble(args[n].substring(6));
		}catch(NumberFormatException error){
		    romb=5;
		}
	    }
	    else if(args[n].startsWith("-seed=")){
		try{
		    seed=Long.parseLong(args[n].substring(6));
		}catch(NumberFormatException error){
		    seed=0;
		}
		SEED=true;
	    }
	    else if(args[n].startsWith("-medi=")){
		med[0]=args[n].substring(6).replace('-', ' ');
	    }
	    else if(args[n].startsWith("-mediadef=")){
		mediadef=args[n].substring(10);
	    }
	    else if(args[n].startsWith("-tdir=")){
		tdir=args[n].substring(6)+"/";
	    }
	    else{
		bnum++;
		pbad+=(bnum>1?",":"")+" \""+args[n]+"\"";
		pflag=false;
	    }
	    if(pflag) param+=" "+args[n];
	}

	Output.err.println(Output.version);
	Output.err.println("Running \""+param+"\"");
	if(bnum>0){
	    param+=" (not used:"+pbad+")";
	    pbad=bnum==1?pbad+" is":"s"+pbad+" are";
	    Output.err.println("Warning: Parameter"+pbad+" not recognized");
	}

	if(gdet!=0 && gdet!=1 && gdet!=2){
	    Output.err.println("Warning: gdet is not a valid number");
	    gdet=0;
	}
	if(length<=0){
	    Output.err.println("Warning: length is not a valid number");
	    length=LENGTH;
	}
	if(radius<=0){
	    Output.err.println("Warning: radius is not a valid number");
	    radius=RADIUS;
	}
	if(width<=0){
	    Output.err.println("Warning: width is not a valid number");
	    width=WIDTH;
	}
	if(height<=0){
	    Output.err.println("Warning: height is not a valid number");
	    height=HEIGHT;
	}
	if(!SURF || surf==0) FACE=false;
	if(bspar<1 || bspar>4){
	    Output.err.println("Warning: bs is not a valid number");
	    bspar=1;
	}
	if(pncrs<1 || pncrs>4 || (pncrs==2 && !(muta.equals("mu") || muta.equals("tau")))){
	    Output.err.println("Warning: ph is not a valid number");
	    pncrs=1;
	}
	if(((pncrs==1 || pncrs==2) && (pncbb<1 || pncbb>4)) ||
	    (pncrs==3 && (pncbb<1 || pncbb>2)) || (pncrs==4 && pncbb!=1)){
	    Output.err.println("Warning: bb is not a valid number");
	    pncbb=1;
	}
	if(((pncrs==1 || pncrs==2) && (pncsh!=1)) || ((pncrs>2) && (pncsh<1 || pncsh>2))){
	    Output.err.println("Warning: sh is not a valid number");
	    pncsh=1;
	}
	if(romb<2 || romb>6){
	    Output.err.println("Warning: romb is not a valid number");
	    romb=5;
	}

	if(mediadef!=null){
	    try{
		if(Output.exists(mediadef)){
		    Reader ifile = Output.reader(mediadef);
		    LineNumberReader inf = new LineNumberReader(ifile);
		    int mediacur;
		    String buf;
		    mediaopts="HI  Attaching the media definition file:\n";
		    while((buf=inf.readLine())!=null){
			mediaopts+="HI  ! "+buf+"\n";
			StringTokenizer st = new StringTokenizer(buf);
			if(!st.hasMoreTokens()) continue;
			String taux=st.nextToken();
			if(taux.charAt(0)=='#') continue;
			if(medianum==mediamax){
			    Output.err.println("Warning: Number of defined media will not exceed "+mediamax);
			    break;
			}
			boolean sflag=false;
			if("all".equals(taux.toLowerCase())){
			    mediacur=0;
			    medt[mediacur]=0;
			}
			else if("sphere".equals(taux.toLowerCase())){
			    mediacur=medianum;
			    medt[mediacur]=1;
			    sphz[mediacur]=Double.parseDouble(st.nextToken());
			    sphr[mediacur]=Double.parseDouble(st.nextToken()); if(sphr[mediacur]<0){ sphr[mediacur]*=-1; sflag=true; }
			    if(sflag) medt[mediacur]*=-1;
			    for(i=0; i<medianum; i++) if(medt[i]==medt[mediacur])
				if(sphz[i]==sphz[mediacur] && sphr[i]==sphr[mediacur]) break;
			    if(i==medianum){ medianum++; mediaseg+=sflag?4:2; }
			}
			else if("box".equals(taux.toLowerCase())){
			    mediacur=medianum;
			    medt[mediacur]=2;
			    boxx[mediacur]=Double.parseDouble(st.nextToken());
			    boxy[mediacur]=Double.parseDouble(st.nextToken());
			    boxz[mediacur]=Double.parseDouble(st.nextToken());
			    boxl[mediacur]=Double.parseDouble(st.nextToken()); if(boxl[mediacur]<0){ boxl[mediacur]*=-1; sflag=true; }
			    boxw[mediacur]=Double.parseDouble(st.nextToken()); if(boxw[mediacur]<0){ boxw[mediacur]*=-1; sflag=true; }
			    boxh[mediacur]=Double.parseDouble(st.nextToken()); if(boxh[mediacur]<0){ boxh[mediacur]*=-1; sflag=true; }
			    if(sflag) medt[mediacur]*=-1;
			    for(i=0; i<medianum; i++) if(medt[i]==medt[mediacur])
				if(boxx[i]==boxx[mediacur] && boxy[i]==boxy[mediacur] && boxz[i]==boxz[mediacur] &&
				   boxl[i]==boxl[mediacur] && boxw[i]==boxw[mediacur] && boxh[i]==boxh[mediacur]) break;
			    if(i==medianum){ medianum++; mediaseg+=sflag?4:2; }
			}
			else if("cyl".equals(taux.toLowerCase())){
			    mediacur=medianum;
			    medt[mediacur]=3;
			    cylz[mediacur]=Double.parseDouble(st.nextToken());
			    cylr[mediacur]=Double.parseDouble(st.nextToken()); if(cylr[mediacur]<0){ cylr[mediacur]*=-1; sflag=true; }
			    cyll[mediacur]=Double.parseDouble(st.nextToken()); if(cyll[mediacur]<0){ cyll[mediacur]*=-1; sflag=true; }
			    if(sflag) medt[mediacur]*=-1;
			    for(i=0; i<medianum; i++) if(medt[i]==medt[mediacur])
				if(cylz[i]==cylz[mediacur] && cylr[i]==cylr[mediacur] && cyll[i]==cyll[mediacur]) break;
			    if(i==medianum){ medianum++; mediaseg+=sflag?4:2; }
			}
			else if("plane".equals(taux.toLowerCase())){
			    mediacur=medianum;
			    medt[mediacur]=4;
			    boxx[mediacur]=Double.parseDouble(st.nextToken());
			    boxy[mediacur]=Double.parseDouble(st.nextToken());
			    boxz[mediacur]=Double.parseDouble(st.nextToken());
			    boxl[mediacur]=Double.parseDouble(st.nextToken());
			    boxw[mediacur]=Double.parseDouble(st.nextToken());
			    boxh[mediacur]=Double.parseDouble(st.nextToken());
			    for(i=0; i<medianum; i++) if(medt[i]==4)
				if(boxx[i]==boxx[mediacur] && boxy[i]==boxy[mediacur] && boxz[i]==boxz[mediacur] &&
				   boxl[i]==boxl[mediacur] && boxw[i]==boxw[mediacur] && boxh[i]==boxh[mediacur]) break;
			    if(i==medianum){ medianum++; mediaseg+=2; }
			}
			else continue;
			medn[mediacur]=Integer.parseInt(st.nextToken())>0?true:false;
			vcut[0][mediacur]=Double.parseDouble(st.nextToken());
			ecut[0][mediacur]=Double.parseDouble(st.nextToken());
			conti[0][mediacur]=Integer.parseInt(st.nextToken())>0?true:false;
			vcut[1][mediacur]=Double.parseDouble(st.nextToken());
			ecut[1][mediacur]=Double.parseDouble(st.nextToken());
			conti[1][mediacur]=Integer.parseInt(st.nextToken())>0?true:false;
			rho[mediacur]=Double.parseDouble(st.nextToken());
			taux=st.nextToken();
			while(st.hasMoreTokens()) taux=taux+" "+st.nextToken();
			med[mediacur]=taux;
		    }
		    ifile.close();
		}
		else Output.err.println("Error: File "+mediadef+" cannot be found");
	    }catch(Exception error){
		Output.err.println("Error: File "+mediadef+" is corrupt: "+error.toString());
		medianum=1; mediaseg=1; medn[0]=true;
	    }
	    mediaopts+="HI  total number of defined media is "+medianum+"\n";
	}

	fnm = new int[mediaseg];
	fnx = new double[mediaseg];

	options=(gdet==2?"HI  radius = "+Output.f(radius):"HI  length = "+Output.f(length)+(gdet==1?" m  width = "+Output.f(width)+" m  height = "+Output.f(height):" m  radius = "+Output.f(radius)))+" m  medium = \""+med[0]+"\"\n";
	if(SURF){
	    options+="HI  propagating to surface z = "+Output.f(surf)+" m";
	    if(FACE){ if(surf<0) options+=" (up only)"; else if(surf>0) options+=" (down only)"; }
	    options+="\n";
	}
	for(i=0; i<medianum; i++){
	    options+="HI  vcut = "+Output.f(vcut[0][i])+"/"+Output.f(vcut[1][i])+"  ecut = "+Output.f(ecut[0][i])+"/"+Output.f(ecut[1][i])+" MeV";
	    if(i==0){
		if(SDEC) options+="  sdec";
		if(RECC) options+="  recc";
		if(timef) options+="  time";
	    }
	    if(conti[0][i]) options+="  cont=1";
	    if(conti[1][i]) options+="  cont=2";
	    if(i==0){
		options+="  bs="+bspar;
		options+="  ph="+pncrs;
		options+="  bb="+pncbb;
		options+="  sh="+pncsh;
		if(lpmef) options+="  lpm";
		if(scatt) options+="  scat";
		options+="  romb="+romb;
	    }
	    else options+="  medium = \""+med[i]+"\"";
	    if(rho[i]!=1) options+="  rho="+Output.f(rho[i]);
	    options+="\n";
	}
	options+=mediaopts;
	if(rfix) options+="HI  Using smart density factor handling, some loss of precision may occur\n";

	Output.AMASIM=amasim;
	PhysicsModel.elow=elow;
	PhysicsModel.ebig=ebig;
	Propagate.g=romb;
	if(SEED){
	    options+="HI  random number generator seed is set at "+seed+"\n";
	    Propagate.rand.setSeed(seed);
	}
	Output.raw=raw;
	Output.err.print(options);

	p1 = new Propagate[medianum];
	p2 = new Propagate[medianum];
	p3 = new Propagate[medianum];

	String settings[] = new String[medianum];

	for(i=0; i<medianum; i++){
	    settings[i]=Output.f(vcut[0][i])+" "+Output.f(ecut[0][i])+" "+conti[0][i]+" "+Output.f(vcut[1][i])+" "+Output.f(ecut[1][i])+" "+conti[1][i];
	    settings[i]+=(frho?" ":(" "+Output.f(rho[i])+" "))+med[i].toLowerCase();
	    for(j=0; j<i; j++) if(settings[j].equals(settings[i])) break;
	    if(i==j){
		p1[i] = new Propagate(med[i], ecut[0][i], vcut[0][i], muta, frho?1:rho[i]);
		p2[i] = new Propagate(med[i], ecut[1][i], vcut[1][i], muta, frho?1:rho[i]);
		p3[i] = new Propagate(med[i], -1.0, -1.0, muta, frho?1:rho[i]);

		if(Output.RecDec) p1[i].sdec=SDEC;  // To enable stopped muon decay
		if(Output.RecDec) p3[i].sdec=SDEC;  // To enable stopped muon decay
		p1[i].contiCorr=conti[0][i];  // To randomize the continuous energy losses, use only for small vcut
		p1[i].exactTime=timef;    // To compute local time of the particle exactly
		p1[i].s.lpm=lpmef;        // Enable lpm and dielectric suppression effects
		p1[i].s.b.form=bspar;     // Choose parametrization of the bremsstrahlung cross section
		p1[i].s.n.form=pncrs;     // Choose parametrization of the photon-nucleon cross section
		p1[i].s.n.bb=pncbb;       // Choose parametrization of the photon-nucleon cross section
		p1[i].s.n.shadow=pncsh;   // Choose parametrization of the photon-nucleon cross section
		p1[i].molieScat=scatt;    // To enable Moliere scattering of the muon
		p1[i].s.ci=crsci;         // Cross section multiplicative modifier
		p1[i].s.cb=crscb;         // Cross section multiplicative modifier
		p1[i].s.cp=crscp;         // Cross section multiplicative modifier
		p1[i].s.ce=crsce;         // Cross section multiplicative modifier
		p1[i].s.cd=crscd;         // Cross section multiplicative modifier
		p2[i].contiCorr=conti[1][i];  // To randomize the continuous energy losses, use only for small vcut
		p2[i].sdec=SDEC;          // To enable stopped muon decay
		p2[i].recc=RECC;          // To enable printout of continuous energy losses
		p2[i].exactTime=timef;    // To compute local time of the particle exactly
		p2[i].s.lpm=lpmef;        // Enable lpm and dielectric suppression effects
		p2[i].s.b.form=bspar;     // Choose parametrization of the bremsstrahlung cross section
		p2[i].s.n.form=pncrs;     // Choose parametrization of the photon-nucleon cross section
		p2[i].s.n.bb=pncbb;       // Choose parametrization of the photon-nucleon cross section
		p2[i].s.n.shadow=pncsh;   // Choose parametrization of the photon-nucleon cross section
		p2[i].molieScat=scatt;    // To enable Moliere scattering of the muon
		p2[i].s.ci=crsci;         // Cross section multiplicative modifier
		p2[i].s.cb=crscb;         // Cross section multiplicative modifier
		p2[i].s.cp=crscp;         // Cross section multiplicative modifier
		p2[i].s.ce=crsce;         // Cross section multiplicative modifier
		p2[i].s.cd=crscd;         // Cross section multiplicative modifier
		p3[i].s.lpm=lpmef;        // Enable lpm and dielectric suppression effects
		p3[i].s.b.form=bspar;     // Choose parametrization of the bremsstrahlung cross section
		p3[i].s.n.form=pncrs;     // Choose parametrization of the photon-nucleon cross section
		p3[i].s.n.bb=pncbb;       // Choose parametrization of the photon-nucleon cross section
		p3[i].s.n.shadow=pncsh;   // Choose parametrization of the photon-nucleon cross section
		p3[i].s.ci=crsci;         // Cross section multiplicative modifier
		p3[i].s.cb=crscb;         // Cross section multiplicative modifier
		p3[i].s.cp=crscp;         // Cross section multiplicative modifier
		p3[i].s.ce=crsce;         // Cross section multiplicative modifier
		p3[i].s.cd=crscd;         // Cross section multiplicative modifier

		p1[i].interpolate("all", tdir+dtct);
		if(medn[i]) p2[i].interpolate("all", tdir+dtct);
		p3[i].interpolate("all", tdir+dtct);
	    }
	    else{
		p1[i]=p1[j];
		p2[i]=p2[j];
		p3[i]=p3[j];
		if(medn[i]) if(!p2[i].jt) p2[i].interpolate("all", tdir+dtct);
	    }
	    if(!frho) rho[i]=1;
	    srho[i]=p1[i].m.Ro*rho[i];
	}
	return true;
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * History lines for the f2k file streams.
     */

    public void history(){
	Output.out.println("HI  "+Output.version+"\nHI  "+param);
	Output.out.print(options);
	Output.out.println("HI  All "+muta+"'s must have energies from "+
			   Output.f(p1[0].p.low*1.e-3)+" GeV to "+Output.f(PhysicsModel.ebig*1.e-3)+" GeV");
	if(USFI) options=" for fixed z = "+Output.f(zset)+" m"; else options="";
	if(USER) Output.out.println("HI  User line "+usna+" will be enabled"+options);
	else{ Output.out.println("HI  Use \"-user\" option to enable "+usna+" user line"); rfix=false; }
	if(rfix) Output.out.println("HI  enforcing compliance with rdmc: exactly one "+usna+" user line per event");
	if(amasim) Output.out.println("HI  some particles will be substituted with AMASIM-compatible types");
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Define user line block for the f2k file streams.
     */

    public void defUSER(){
	if(USER){
	    Output.out.println("USER_DEF "+usna+" NR E_INI E_CPD E_IN E_OUT CDP_X CDP_Y CDP_Z Z_IN Z_OUT");
	    if(rw!=0) Output.out.println("USER_DEF ev_wght WEIGHT DIST_R");
	}
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize user line block for the f2k file streams.
     */

    public void iniUSER(){
	if(USER){ userbf.clear(); imax=0; emax=0; }
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Record user line block for the f2k file streams.
     */

    public void recUSER(){
	userbf.addElement(igen+" "+Output.f(e)+" "+(ec>0?Output.f(ec):"?")+" "+
			  (e1!=0?Output.f(e1):"?")+" "+(e2!=0?Output.f(e2):"?")+" "+
			  Output.f(nx)+" "+Output.f(ny)+" "+Output.f(nz)+" "+
			  (h1>0?Output.f(z1):"?")+" "+(h2>0?Output.f(z2):"?"));
	if(rfix) if(e>emax){ imax=userbf.size()-1; emax=e; }
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Output user line block for the f2k file streams.
     */

    public void outUSER(){
	if(USER){
	    if(rfix){
		if(userbf.size()==0) Output.out.println("US "+usna+" ? ? ? ? ? ? ? ? ? ?");
		else Output.out.println("US "+usna+" "+userbf.elementAt(imax));
	    }
	    else for(int i=0; i<userbf.size(); i++) Output.out.println("US "+usna+" "+userbf.elementAt(i));
	    if(rw!=0) Output.out.println("US ev_whgt "+fw+" "+hw);
	}
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize particle buffers for the f2k file streams.
     */

    public void initBufs(){
	hist = new StringBuffer(Output.HISTSIZE);
	userbf = new Vector(Output.HISTSIZE);
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Main parser for the f2k file streams.
     */

    public void mmcf2k(String[] args){
	Output.I3flag=false; dtct=".amanda";
	if(!setp(args)) return;

	Output.err.println(" ---  *** Enter your input in F2000 format now ***  --- ");

	try{
	    LineNumberReader file = new LineNumberReader(new InputStreamReader(Output.in));
	    StringTokenizer st;
	    String line, taux;
	    Vector buffer = new Vector(Output.HISTSIZE);
	    initBufs();
	    boolean tbegin=false;

	    int i, ipar;

	    if((line=file.readLine())!=null){
		Output.out.println(line);
		history();
	    }

	    while((line=file.readLine())!=null){
		st = new StringTokenizer(line);
		if(!st.hasMoreTokens()) continue;
		taux=st.nextToken();
		if("TBEGIN".equals(taux)){
		    defUSER();
		    Output.out.println(line);
		    tbegin=true;
		    continue;
		}
		if(!tbegin){
		    Output.out.println(line);
		    continue;
		}
		buffer.addElement(line);
		if("TR".equals(taux)){
		    igen=Integer.parseInt(st.nextToken());
		    if(igen>gens) gens=igen;
		    prnt=st.nextToken();
		    try{
			ipar=Integer.parseInt(prnt);
			if(ipar>gens) gens=ipar;
		    }catch(NumberFormatException error){
			ipar=-1;
		    }
		}
		else if("EE".equals(taux) || "END".equals(taux)){
		    iniUSER();
		    if(rw!=0) dw=true;

		    for(i=0; i<buffer.size(); i++){
			line=(String)buffer.elementAt(i);
			st = new StringTokenizer(line);
			taux=st.nextToken();
			if("TR".equals(taux)){

			    igen=Integer.parseInt(st.nextToken());
			    prnt=st.nextToken();
			    type=st.nextToken();
			    if(muta.equals(type) || (muta+"-").equals(type) || (muta+"+").equals(type)){

				double x, y, z, th, phi, l, e, t;

				x=Double.parseDouble(st.nextToken());
				y=Double.parseDouble(st.nextToken());
				z=Double.parseDouble(st.nextToken());
				th=Double.parseDouble(st.nextToken());
				phi=Double.parseDouble(st.nextToken());
				if(lfix){
				    st.nextToken();
				    l=0;
				}
				else try{
				    l=Double.parseDouble(st.nextToken());
				    if(l<0) l=0;
				    else l*=1.e2;
				}catch(NumberFormatException error){
				    l=0;
				}
				e=Double.parseDouble(st.nextToken());
				t=Double.parseDouble(st.nextToken());

				prop(x, y, z, th, phi, l, e, t);

			    }
			    else{
				type="TR "+igen+" "+prnt+" "+type+st.nextToken("\n");
				Output.out.println(type);
			    }
			}
			else if("EM".equals(taux)){
			    events++;
			    Output.out.println(line);
			}
			else if("EE".equals(taux)){
			    outUSER();
			    Output.out.println(line);
			}
			else Output.out.println(line);
		    }

		    gens=0;
		    buffer.clear();
		    if(rw!=0) dw=false;
		}

	    }

	    if(events!=0) vertices/=events;
	    Output.err.println("events read "+events+" tracks looped "+tracks+
			       " average vertices "+vertices+" missed volume "+missed);
	}catch(Exception error){
	    Output.err.println("Program finished with exception: "+error.toString());
	    throw new mmcException("input error");
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Main propagator routine.
     */

    public double prop(int igen, int gens, double x, double y, double z, double th, double phi, double l, double e, double t){
	this.igen=igen;
	this.gens=gens;
	return prop(x, y, z, th, phi, l, e, t);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Main propagator routine.
     */

    public double prop(double x, double y, double z, double th, double phi, double l, double e, double t){
	this.e=e;
	tracks++;
	int j, m, mi, mf, md, mj;

	Propagate pp=null;
	double r1=0, r2=0, r3, rr, h3, n2=0, aux;
	double dx, dy, dz;
	double result;

	double sinth, sinph, costh, cosph;

	sinph=Math.sin(phi*Math.PI/180);
	cosph=Math.cos(phi*Math.PI/180);
	sinth=Math.sin(th*Math.PI/180);
	costh=Math.cos(th*Math.PI/180);

	dx=-sinth*cosph;
	dy=-sinth*sinph;
	dz=-costh;

	switch(gdet){
	case 1:
	    setbox(x, y, z, dx, dy, dz, length, width, height);
	    break;
	case 2:
	    setsph(x, y, z-radius, cosph, sinph, costh, sinth, radius);
	    break;
	case 0:
	default:
	    setcyl(x, y, z, cosph, sinph, costh, sinth, length, radius);
	}
	h1=hx1; h2=hx2;

	if(h1<0) h1=0;
	if(h2<0) h2=0;
	if(h1==h2) missed++;
	h3=BIG;
	if(dw){ hw=h1+(h2-h1)*Propagate.rand.nextDouble(); fw=1; }

	if(SURF){
	    if(costh!=0) aux=(z-surf)/costh;
	    else aux=BIG;
	    if(aux<0) aux=BIG;
	    if(FACE) if((z-surf)*surf<0) aux=0;
	    if(h1>aux) h1=aux;
	    if(h2>aux) h2=aux;
	    if(h3>aux) h3=aux;
	}

	if(USER){
	    if(USFI){
		if(costh!=0) n2=(z-zset)/costh;
		else n2=0;
	    }
	    else n2=x*sinth*cosph+y*sinth*sinph+z*costh;
	    if(n2<h1) flag=1;
	    else if(n2<h2) flag=2;
	    else if(n2<h3) flag=3;
	    else flag=0;
	}

	fnu=1;
	fnm[0]=0;
	fnx[0]=BIG;
	for(m=1; m<medianum; m++){
	    switch(Math.abs(medt[m])){
	    case 1:
		setsph(x, y, z-sphz[m], cosph, sinph, costh, sinth, sphr[m]);
		break;
	    case 2:
		setbox(x-boxx[m], y-boxy[m], z-boxz[m], dx, dy, dz, boxl[m], boxw[m], boxh[m]);
		break;
	    case 3:
		setcyl(x, y, z-cylz[m], cosph, sinph, costh, sinth, cyll[m], cylr[m]);
		break;
	    case 4:
		setpln(x-boxx[m], y-boxy[m], z-boxz[m], dx, dy, dz, boxl[m], boxw[m], boxh[m]);
		break;
	    default: hx1=0; hx2=0;
	    }
	    for(int ij=0; ij<(medt[m]>0?1:2); ij++){
		if(medt[m]>0){ r1=hx1; r2=hx2; }
		else{
		    if(ij==0){ r1=0; r2=hx1; }
		    else{ r1=hx2; r2=BIG; }
		}
		if(r2>r1){
		    md=0;
		    mi=0; if(r1>0){ md++; for(mi=0; mi<fnu; mi++) if(fnx[mi]>r1) break; }
		    mf=0; if(r2>0){ md++; for(mf=mi; mf<fnu; mf++) if(fnx[mf]>r2) break; }
		    md+=mi-mf;
		    if(md>0) for(mj=fnu-1; mj>=mf; mj--){
			fnx[mj+md]=fnx[mj];
			fnm[mj+md]=fnm[mj];
		    }
		    else if(md<0) for(mj=mf; mj<fnu; mj++){
			fnx[mj+md]=fnx[mj];
			fnm[mj+md]=fnm[mj];
		    }
		    fnu+=md;
		    if(r1>0){
			fnx[mi]=r1;
			mi++;
		    }
		    if(r2>0){
			fnx[mi]=r2;
			fnm[mi]=m;
		    }
		}
	    }
	}
	if(DEBUG){
	    for(j=0; j<fnu; j++) Output.out.println(j+"\t"+Output.f(fnx[j])+"\t"+fnm[j]+"\t"+p1[fnm[j]].m.name+"\t"+Output.f(rho[fnm[j]]));
	    Output.out.println();
	}

	if(Output.I3flag) I3hist.clear(); else hist.setLength(0);
	if(USER){ e1=0; e2=0; ec=0; elost=0; }
	pp=null; result=e; rr=0;

	for(j=0; j<fnu; j++){
	    m=fnm[j]; aux=fnx[j];
	    r1=Math.min(h1, aux);
	    r2=Math.min(h2, aux);
	    r3=Math.min(h3, aux);
	    p1[m].rho=rho[m];
	    p2[m].rho=rho[m];
	    p3[m].rho=rho[m];
	    if(pp==null) result=propagateTo(r1, e, n2, 1, p1[m], igen, gens,
					    t*1.e-9, x*1.e2, y*1.e2, z*1.e2, 180-th, phi<180?phi+180:phi-180);
	    else result=propagateTo(r1-rr, result, n2-rr, 1, p1[m], 0, 0,
				    pp.p.t, pp.p.x, pp.p.y, pp.p.z, pp.p.theta, pp.p.phi);
	    pp=p1[m]; l+=pp.p.r; if(USER) if(r1>=rr && r1==h1){ e1=result; if(Output.I3flag) pI.ti=pp.p.t; }
	    if(Output.RecDec){ if(Output.I3flag) I3hist.addAll(pp.o.I3hist); else hist.append(pp.o.history); gens=pp.o.gens; }
	    if(result>0){
		if(r1>rr) rr=r1; if(USER) elost+=result;
		if(dw){ p2[m].dw=true; p2[m].rw=rw; p2[m].hw=(hw-rr)*1.e2; }
		result=propagateTo(r2-rr, result, n2-rr, 2, p2[m], igen, gens,
				   pp.p.t, pp.p.x, pp.p.y, pp.p.z, pp.p.theta, pp.p.phi);
		pp=p2[m]; l+=pp.p.r; if(USER) if(r2>=rr && r2==h2){ e2=result; if(Output.I3flag) pI.tf=pp.p.t; }
		if(dw){ if(!p2[m].dw){ dw=false; fw=p2[m].rw; hw=p2[m].hw*1.e-2+rr; } else p2[m].dw=false; p2[m].rw=0; p2[m].hw=0; }
		if(Output.I3flag) I3hist.addAll(pp.o.I3hist); else hist.append(pp.o.history);
		vertices+=pp.o.gens-gens; gens=pp.o.gens;
		if(result>0){
		    if(r2>rr) rr=r2; if(USER){ elost-=result; if(elost<0) elost=0; }
		    result=propagateTo(r3-rr, result, n2-rr, 3, p3[m], igen, gens,
				       pp.p.t, pp.p.x, pp.p.y, pp.p.z, pp.p.theta, pp.p.phi);
		    pp=p3[m]; l+=pp.p.r;
		    if(Output.RecDec){ if(Output.I3flag) I3hist.addAll(pp.o.I3hist); else hist.append(pp.o.history); gens=pp.o.gens; }
		    if(result>0) if(r3>rr) rr=r3;
		}
	    }
	    if(result<=0) break;
	}

	if(USER){
	    if(h1==0) e1=0;
	    if(h2==0) e2=0;
	    if(e1<0) e1-=rr;
	    if(e2<0) e2-=rr;
	    z1=z-costh*h1;
	    z2=z-costh*h2;
	    nx=x-sinth*cosph*n2;
	    ny=y-sinth*sinph*n2;
	    nz=z-costh*n2;
	    if(!Output.I3flag) recUSER();
	    else{
		pI.xi=x-sinth*cosph*h1;
		pI.yi=y-sinth*sinph*h1;
		pI.zi=z-costh*h1;
		pI.xf=x-sinth*cosph*h2;
		pI.yf=y-sinth*sinph*h2;
		pI.zf=z-costh*h2;
		pI.xc=nx;
		pI.yc=ny;
		pI.zc=nz;
		pI.Ei=e1;
		pI.Ef=e2;
		pI.Ec=ec>0?ec:0;
		pI.Elost=elost;
	    }
	}

	if(SURF){
	    x=pp.p.x*1.e-2;
	    y=pp.p.y*1.e-2;
	    z=pp.p.z*1.e-2;
	    th=180-pp.p.theta;
	    phi=pp.p.phi; phi=phi<180?phi+180:phi-180;
	    if(result>0){ e=result; l=-BIG; }
	    else{ e=0; l=0; }
	    t=pp.p.t*1.e9;
	}
	if(!Output.I3flag){
	    Output.out.println("TR "+igen+" "+prnt+" "+type+" "+Output.f(x)+" "+Output.f(y)+" "+Output.f(z)+" "+
			       Output.f(th)+" "+Output.f(phi)+" "+(l>=0?Output.f(l*1.e-2):"?")+" "+
			       Output.f(e)+" "+Output.f(t));
	    Output.out.print(hist);
	}

	return l;
    }


    //----------------------------------------------------------------------------------------------------//

    private double hx1, hx2;

    /**
     * Calculate points of intersection with the cylinder.
     */

    private void setcyl(double x, double y, double z, double cosph, double sinph, double costh,double sinth, double length, double radius){

	double aux, r1, r2;
	double b=x*cosph+y*sinph;
	double d=b*b+radius*radius-x*x-y*y;
	if(d>0){
	    d=Math.sqrt(d);
	    if(costh!=0){
		hx1=(z-length/2)/costh;
		hx2=(z+length/2)/costh;
		if(hx1>hx2){ aux=hx1; hx1=hx2; hx2=aux; }
	    }
	    if(sinth!=0){
		r1=(b-d)/sinth;
		r2=(b+d)/sinth;
		if(r1>r2){ aux=r1; r1=r2; r2=aux; }
		if(costh==0){
		    if(z>-length/2 && z<length/2){ hx1=r1; hx2=r2; }
		    else{ hx1=0; hx2=0; }
		}
		else{
		    if(hx1>=r2 || hx2<=r1){ hx1=0; hx2=0; }
		    else{ hx1=Math.max(r1, hx1); hx2=Math.min(r2, hx2); }
		}
	    }
	}
	else{ hx1=0; hx2=0; }
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculate points of intersection with the box.
     */

    private void setbox(double x, double y, double z, double dx, double dy, double dz, double length, double width, double height){
	boolean top=false, bottom=false, tbflag=false;
	double hx;
	if(dx!=0){
	    if(!tbflag){
		hx=-(x-length/2)/dx;
		if(Math.abs(y+hx*dy)<=width/2 && Math.abs(z+hx*dz)<=height/2){
		    if(!top){ hx1=hx; top=true; }
		    else if(hx!=hx1 && !bottom){ hx2=hx; bottom=true; }
		    else tbflag=true;
		}
	    }
	    if(!tbflag){
		hx=-(x+length/2)/dx;
		if(Math.abs(y+hx*dy)<=width/2 && Math.abs(z+hx*dz)<=height/2){
		    if(!top){ hx1=hx; top=true; }
		    else if(hx!=hx1 && !bottom){ hx2=hx; bottom=true; }
		    else tbflag=true;
		}
	    }
	}
	if(dy!=0){
	    if(!tbflag){
		hx=-(y-width/2)/dy;
		if(Math.abs(x+hx*dx)<=length/2 && Math.abs(z+hx*dz)<=height/2){
		    if(!top){ hx1=hx; top=true; }
		    else if(hx!=hx1 && !bottom){ hx2=hx; bottom=true; }
		    else tbflag=true;
		    }
	    }
	    if(!tbflag){
		hx=-(y+width/2)/dy;
		if(Math.abs(x+hx*dx)<=length/2 && Math.abs(z+hx*dz)<=height/2){
		    if(!top){ hx1=hx; top=true; }
		    else if(hx!=hx1 && !bottom){ hx2=hx; bottom=true; }
		    else tbflag=true;
		}
	    }
	}
	if(dz!=0){
	    if(!tbflag){
		hx=-(z-height/2)/dz;
		if(Math.abs(x+hx*dx)<=length/2 && Math.abs(y+hx*dy)<=width/2){
		    if(!top){ hx1=hx; top=true; }
		    else if(hx!=hx1 && !bottom){ hx2=hx; bottom=true; }
		    else tbflag=true;
		}
	    }
	    if(!tbflag){
		hx=-(z+height/2)/dz;
		if(Math.abs(x+hx*dx)<=length/2 && Math.abs(y+hx*dy)<=width/2){
		    if(!top){ hx1=hx; top=true; }
		    else if(hx!=hx1 && !bottom){ hx2=hx; bottom=true; }
		    else tbflag=true;
		}
	    }
	}
	if(!top) hx1=0;
	if(!bottom) hx2=0;
	if(hx1>hx2){ hx=hx1; hx1=hx2; hx2=hx; }
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculate points of intersection with the sphere.
     */

    private void setsph(double x, double y, double z, double cosph, double sinph, double costh, double sinth, double radius){
	double b=(x*cosph+y*sinph)*sinth+(z+radius)*costh;
	double d=b*b-(x*x+y*y+z*(z+2*radius));
	if(d>0){
	    d=Math.sqrt(d);
	    hx1=b-d;
	    hx2=b+d;
	}
	else{
	    hx1=0;
	    hx2=0;
	}
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculate points of intersection with the box.
     */

    private void setpln(double x, double y, double z, double dx, double dy, double dz, double nx, double ny, double nz){
	double b=-(x*nx+y*ny+z*nz);
	double d=dx*nx+dy*ny+dz*nz;
	if(b>=0 || b*d>0){
	    if(b>=0){
		hx1=0;
		hx2=d!=0?b/d:BIG;
	    }
	    else{
		hx1=d!=0?b/d:0;
		hx2=BIG;
	    }
	}
	else{
	    hx1=0;
	    hx2=0;
	}
    }


    //----------------------------------------------------------------------------------------------------//

    public boolean dw=false;
    public double rw=0, hw=0, fw=1;

    private int flag=0;
    private double e;

    /**
     * user-block variable: z-coordinate of entry point into the detector cylinder (Z_IN)
     */

    public double z1=0;

    /**
     * user-block variable: z-coordinate of exit point from the detector cylinder (Z_OUT)
     */

    public double z2=0;

    /**
     * user-block variable: distance to the point of entry into the detector cylinder (R_IN)
     */

    public double h1=0;

    /**
     * user-block variable: distance to the point of exit from the detector cylinder (R_OUT)
     */

    public double h2=0;

    /**
     * user-block variable: x-coordinate of closest approach point (CPD_X) or -user=z point
     */

    public double nx=0;

    /**
     * user-block variable: y-coordinate of closest approach point (CPD_Y) or -user=z point
     */

    public double ny=0;

    /**
     * user-block variable: z-coordinate of closest approach point (CPD_Z) or -user=z point
     */

    public double nz=0;

    /**
     * user-block variable: Energy at the entry point into the detector cylinder if positive (E_IN)
     */

    public double e1=0;

    /**
     * user-block variable: Energy at the exit point from the detector cylinder if positive (E_OUT)
     */

    public double e2=0;

    /**
     * user-block variable: Energy at the CPD point if positive (E_CPD)
     */

    public double ec=0;

    /**
     * Energy lost inside the detector cylinder
     */

    public double elost=0;

    //----------------------------------------------------------------------------------------------------//

    /**
     * checks the flag, corrects the units, initializes the particle and calls the propagator
     */

    private double propagateTo(double h, double e, double n, int i, Propagate p, int igen, int gens,
				      double time, double x, double y, double z, double theta, double phi){
	double result;
	switch(i){
	case 1: p.o.initDefault(igen, gens, type, time, x, y, z, theta, phi); break;
	case 2: p.o.initF2000(igen, gens, type, time, x, y, z, theta, phi); break;
	case 3: p.o.initDefault(igen, gens, type, time, x, y, z, theta, phi); break;
	}
	if(flag==i && n>=0 && n<h){
	    result=p.propagateTo(h*1.e2, e*1.e3, n*1.e2);
	    ec=p.getPropEc();
	    ec=ec>0?ec*1.e-3:ec*1.e-2;
	    if(Output.I3flag) pI.tc=p.getPropTc();
	}
	else result=p.propagateTo(h*1.e2, e*1.e3);
	return result>0?result*1.e-3:result*1.e-2;
    }

}
