package gen;
import mmc.*;
import tfa.*;
import java.util.Random;
import java.util.Vector;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.util.StringTokenizer;

/**
 * This is the main lepton generator class.
 */

public class AtmFlux extends PhysicsModel{

    public double D=1730;
    public double length, radius;
    private double R=400;
    private double tf[], af[], sF[], tF=0;
    private double nF=1, mF=0, xn=1;
    private double tT=1;
    private double t0=0;
    private long evt=0;
    private boolean OSC=true;
    private boolean GEN=false;
    private boolean PROP=false;

    private double R0=EarthModel.R0*1.e3;

    private int gens, lept=0;
    private double px, py, pz, th, phi;
    public double elost, eini, thini;
    public String Name;

    public EarthModel eM;
    private IntFlux iF[];
    private NeutrinoTot nT;
    private Amanda mP, tP;
    private Integral I;
    private Random rand;

    private Vector I3hist;
    private boolean f2k=true;
    private double Ecut[], A[], G[], a[], g[];
    private String tdir="";
    private String dtct="";
    private String name[] = {"mu-", "mu+", "nu_mu", "~nu_mu", "nu_e", "~nu_e", "nu_tau", "~nu_tau", "e-", "e+", "tau-", "tau+", "hadr"};


    //----------------------------------------------------------------------------------------------------//

    /**
     * Class initializer, command-line option parser. Call with "-help" to list all options.
     */

    public long initAtmFlux(String[] args, boolean f2k){
	this.f2k=f2k;
	String param="AtmFlux";
	int i;
	long num=10;

	double Emu=300;
	double Enm=100;
	double Ene=100;
	double Ent=100;
	double Amu=1, Anm=1, Ane=1;
	double Gmu=0, Gnm=0, Gne=0;
	double ant=0, anm=0, ane=0;
	double gnt=2, gnm=2, gne=2;
	double prompt=0;

	int w=1, f=0, m=-1;
	double h0=-114.8, z0=2.834, b0=0.024;
	int sM=0;
	double Efr=10;
	double gE=0, gD=0.2;

	long seed=0;
	boolean SEED=false;

	int bnum=0;
	String pbad="";
	boolean pflag;

	for(i=0; i<args.length; i++){
	    pflag=true;
	    if(args[i].equals("-help") || args[i].equals("-h") || args[i].equals("--help")){
		Output.out.println("\n"+
"This program generates atmospheric lepton fluxes at 0.1-4000 TeV\n"+
"                       -Emu=[100-1.e5 GeV] muon low energy threshold\n"+
"                       -Enm=[GeV] muon neutrino low energy threshold\n"+
"                       -Ene=[GeV] electron neutrino energy threshold\n"+
"                       -Ent=[GeV] tau neutrino  low energy threshold\n"+
"                       -A[mu/nm/ne]=[normalization  correction]\n"+
"                       -G[mu/nm/ne]=[spectral index correction]\n"+
"                       -a[nm/ne/nt]=[    for the ad-hoc component    ]\n"+
"                       -g[nm/ne/nt]=[1e-7/(GeV sr s cm^2) a(E/GeV)^-g]\n"+
"                       -prompt=[R_c]  ratio of prompt to pion muons\n"+
"                       -R=[radius] of the detector in meters\n"+
"                       -D=[depth]  of the detector in meters\n"+
"                       -N=[number] of events to generate\n"+
"                       -f=[0-1]  method of cos* correction\n"+
"                       -m=[-1/0-3] choose atmosphere model\n"+
"                       -z0=[ground  elevation in km]\n"+
"                       -b0=[bedrock elevation in km]\n"+
"                       -h0=[average production height in km\n"+
"                                   or in g/cm^2 if negative]\n"+
"                       -nF=[neutrino cross sections multiplic. factor]\n"+
"                       -mF=[max. mass overburden/neutrino int. length]\n"+
"                       -sM=[0-4]  Solar modulation:  none/min/mid/max\n"+
"                       -Efr=[cr primary to neutrino energy fraction]\n"+
"                       -gE=[GeV]  artificial geomagnetic cutoff value\n"+
"                       -gD=[decades]  span of the response function\n"+
"                       -noosc   disable mu < - > tau oscillations\n"+
"                       -nlow=[GeV]   lowest neutrino energy\n"+
"                       -gen   use as generator  only\n"+
"                       -prop  use as propagator only\n"+
"                       -seed=[integer] random number generator seed\n"+
"mmc help lines follow:");
		new Amanda().setp(args);
		return -1;
	    }
	    else if(args[i].equals("-gen")){
		GEN=true;
	    }
	    else if(args[i].equals("-prop")){
		PROP=true;
	    }
	    else if(args[i].equals("-noosc")){
		OSC=false;
	    }
	    else if(args[i].startsWith("-Emu=")){
		try{
		    Emu=Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    Emu=600;
		}
	    }
	    else if(args[i].startsWith("-Enm=")){
		try{
		    Enm=Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    Enm=600;
		}
	    }
	    else if(args[i].startsWith("-Ene=")){
		try{
		    Ene=Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    Ene=600;
		}
	    }
	    else if(args[i].startsWith("-Ent=")){
		try{
		    Ent=Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    Ent=600;
		}
	    }
	    else if(args[i].startsWith("-Amu=")){
		try{
		    Amu=Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    Amu=1;
		}
	    }
	    else if(args[i].startsWith("-Anm=")){
		try{
		    Anm=Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    Anm=1;
		}
	    }
	    else if(args[i].startsWith("-Ane=")){
		try{
		    Ane=Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    Ane=1;
		}
	    }
	    else if(args[i].startsWith("-Gmu=")){
		try{
		    Gmu=Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    Gmu=0;
		}
	    }
	    else if(args[i].startsWith("-Gnm=")){
		try{
		    Gnm=Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    Gnm=0;
		}
	    }
	    else if(args[i].startsWith("-Gne=")){
		try{
		    Gne=Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    Gne=0;
		}
	    }
	    else if(args[i].startsWith("-ant=")){
		try{
		    ant=Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    ant=0;
		}
	    }
	    else if(args[i].startsWith("-anm=")){
		try{
		    anm=Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    anm=0;
		}
	    }
	    else if(args[i].startsWith("-ane=")){
		try{
		    ane=Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    ane=0;
		}
	    }
	    else if(args[i].startsWith("-gnt=")){
		try{
		    gnt=Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    gnt=2;
		}
	    }
	    else if(args[i].startsWith("-gnm=")){
		try{
		    gnm=Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    gnm=2;
		}
	    }
	    else if(args[i].startsWith("-gne=")){
		try{
		    gne=Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    gne=2;
		}
	    }
	    else if(args[i].startsWith("-prompt=")){
		try{
		    prompt=Double.parseDouble(args[i].substring(8));
		}catch(Exception error){
		    prompt=0;
		}
	    }
	    else if(args[i].startsWith("-R=")){
		try{
		    R=Double.parseDouble(args[i].substring(3));
		}catch(Exception error){
		    R=10;
		}
	    }
	    else if(args[i].startsWith("-D=")){
		try{
		    D=Double.parseDouble(args[i].substring(3));
		}catch(Exception error){
		    D=0;
		}
	    }
	    else if(args[i].startsWith("-N=")){
		try{
		    num=(long)Double.parseDouble(args[i].substring(3));
		}catch(Exception error){
		    num=1;
		}
	    }
	    else if(args[i].startsWith("-f=")){
		try{
		    f=(int)Double.parseDouble(args[i].substring(3));
		}catch(Exception error){
		    f=0;
		}
	    }
	    else if(args[i].startsWith("-m=")){
		try{
		    m=(int)Double.parseDouble(args[i].substring(3));
		}catch(Exception error){
		    m=-1;
		}
	    }
	    else if(args[i].startsWith("-z0=")){
		try{
		    z0=Double.parseDouble(args[i].substring(4));
		}catch(Exception error){
		    z0=0;
		}
	    }
	    else if(args[i].startsWith("-b0=")){
		try{
		    b0=Double.parseDouble(args[i].substring(4));
		}catch(Exception error){
		    b0=0;
		}
	    }
	    else if(args[i].startsWith("-h0=")){
		try{
		    h0=Double.parseDouble(args[i].substring(4));
		}catch(Exception error){
		    h0=-114.8;
		}
	    }
	    else if(args[i].startsWith("-nF=")){
		try{
		    nF=Double.parseDouble(args[i].substring(4));
		}catch(Exception error){
		    nF=1;
		}
	    }
	    else if(args[i].startsWith("-mF=")){
		try{
		    mF=Double.parseDouble(args[i].substring(4));
		}catch(Exception error){
		    mF=0;
		}
	    }
	    else if(args[i].startsWith("-sM=")){
		try{
		    sM=(int)Double.parseDouble(args[i].substring(4));
		}catch(Exception error){
		    sM=0;
		}
	    }
	    else if(args[i].startsWith("-Efr=")){
		try{
		    Efr=Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    Efr=10;
		}
	    }
	    else if(args[i].startsWith("-gE=")){
		try{
		    gE=Double.parseDouble(args[i].substring(4));
		}catch(Exception error){
		    gE=0;
		}
	    }
	    else if(args[i].startsWith("-gD=")){
		try{
		    gD=Double.parseDouble(args[i].substring(4));
		}catch(Exception error){
		    gD=0.2;
		}
	    }
	    else if(args[i].startsWith("-nlow=")){
		try{
		    nlow=Double.parseDouble(args[i].substring(6))*1.e3;
		}catch(Exception error){
		    nlow=Me;
		}
	    }
	    else if(args[i].equals("-tau")){
		args[i]="-ignore-tau";
	    }
	    else if(args[i].startsWith("-seed=")){
		try{
		    seed=Long.parseLong(args[i].substring(6));
		}catch(Exception error){
		    seed=0;
		}
		SEED=true;
	    }
	    else if(args[i].startsWith("-tdir=")){
		tdir=args[i].substring(6)+"/";
	    }
	    else{
		bnum++;
		pbad+=(bnum>1?",":"")+" \""+args[i]+"\"";
		pflag=false;
	    }
	    if(pflag) param+=" "+args[i];
	}

	Output.err.println(Output.version);
	Output.err.println("Running \""+param+"\"");
	if(bnum>0){
	    param+=" (not used:"+pbad+")";
	    pbad=bnum==1?pbad+" is":"s"+pbad+" are";
	    Output.err.println("Warning: Parameter"+pbad+" not recognized");
	}

	if(Emu<Mmu*1.e-3) Emu=Mmu*1.e-3;
	if(Enm<Me*1.e-3) Enm=Me*1.e-3;
	if(Ene<Me*1.e-3) Ene=Me*1.e-3;
	if(Ent<Me*1.e-3) Ent=Me*1.e-3;

	if(Amu<0) Amu=0;
	if(Anm<0) Anm=0;
	if(Ane<0) Ane=0;
	if(ant<0) ant=0;
	if(anm<0) anm=0;
	if(ane<0) ane=0;

	if(Gmu<-1) Gmu=0;
	if(Gnm<-1) Gnm=0;
	if(Gne<-1) Gne=0;
	if(gnt<=1) gnt=2;
	if(gnm<=1) gnm=2;
	if(gne<=1) gne=2;

	if(Amu==0 && Anm==0 && Ane==0 && ant==0 && anm==0 && ane==0) Amu=1;
	if(GEN) PROP=false;

	if(h0==0) m=-1;
	if(f<0 || f>1) f=0;
	if(m<-1 || m>3) m=-1;
	Output.err.println("Choosing f="+f+" m="+m+" z0="+Output.f(z0)+" km  b0="+
			   Output.f(b0)+" km  h0="+Output.f(h0)+" "+(h0>0?"km":"g/cm^2"));
	Output.err.println("Detector with radius R="+Output.f(R)+" m is at depth="+Output.f(D)+" m");
	Output.err.println("Energy thresholds: Emu="+Output.f(Emu)+" GeV, Enm="+Output.f(Enm)+" GeV, Ene="+Output.f(Ene)+" GeV, Ent="+Output.f(Ent)+" GeV");

	if(nF<=0) nF=1;
	if(nF>1) Output.err.println("Use nF>1 to speed up the calculation for energy thresholds << 10 TeV");
	else if(nF<1) Output.err.println("Setting nF value to "+Output.f(nF)+" < 1. Your time would probably be much better spent elsewhere");

	if(mF<1) mF=0;
	if(GEN || PROP) mF=0;
	if(mF>0) Output.err.println("Using re-weighting technique. This is an approximation! Only use mF>>1, e.g. mF=10");

	if(sM<0) sM=0;
	if(PROP) sM=0;
	if(Efr<1) Efr=10;
	if(sM>0) Output.err.println("Solar modulation enabled with parameters: sM="+sM+", Efr="+Output.f(Efr)+" GeV");
	if(gE<0) gE=0;
	if(gD<=0) gD=0.2;
	if(gE>0) Output.err.println("Geomagnetic cutoff enabled with parameters: gE="+Output.f(gE)+
				    " GeV, gD="+Output.f(gD)+(sM>0?"":", Efr="+Output.f(Efr)+" GeV"));

	Decay.flag=true;
	Output.RecDec=true;
	if(!GEN){
	    mP = new Amanda();
	    tP = new Amanda();
	    mP.dtct=dtct;
	    tP.dtct=dtct;
	    eM = new EarthModel(m<0?0:m, z0, b0, D*1.e-3);
	    nT = new NeutrinoTot();
	}

	rand = new Random();
	R0+=z0*1.e3;

	if(f2k) Output.out.println("V 2000.1.2");

	String ename="";
	if(!GEN){
	    mP.setup(args);
	    if(f2k) mP.history();
	    String argt[] = new String[args.length+1];
	    for(i=0; i<args.length; i++) argt[i]=args[i];
	    argt[args.length]="-tau";
	    tP.setup(argt);
	    if(f2k) tP.history();
	    switch(mP.gdet){
	    case 1:
		length=mP.height;
		radius=mP.length/2;
		break;
	    case 2:
		length=mP.radius*2;
		radius=mP.radius;
		break;
	    case 0:
	    default:
		length=mP.length;
		radius=mP.radius;
	    }

	    eM.num=mP.medianum;
	    eM.mRa = new double[eM.num];
	    eM.mRo = new double[eM.num];
	    eM.mRa[0]=R0*1.e-3;
	    eM.mRo[0]=mP.srho[0];
	    ename="_"+eM.mRo[0];
	    for(i=1; i<eM.num; i++){
		if(mP.medt[i]==1 && Math.abs((mP.sphr[i]-mP.sphz[i])-(R0-D))<R0*halfPrecision){
		    eM.mRa[i]=(mP.sphz[i]+R0-D)*1.e-3;
		    eM.mRo[i]=mP.srho[i];
		    ename+="_"+Output.f(mP.sphz[i])+"-"+Output.f(eM.mRo[i]);
		}
		else eM.mRo[i]=0;
	    }
	}

	String name=tdir+dtct;

	if(!PROP){
	    I = new Integral(iromb, imaxs, iprec);
	    int type[] = {2, 3, 5, 6, 8, 9};
	    tf = new double[6];
	    af = new double[6];
	    sF = new double[6];
	    Ecut = new double[8];
	    A = new double[6];
	    G = new double[6];
	    a = new double[6];
	    g = new double[6];
	    Ecut[0]=Emu; A[0]=Amu; G[0]=Gmu; a[0]=ant*1.e-7; g[0]=gnt;
	    Ecut[1]=Emu; A[1]=Amu; G[1]=Gmu; a[1]=ant*1.e-7; g[1]=gnt;
	    Ecut[2]=Enm; A[2]=Anm; G[2]=Gnm; a[2]=anm*1.e-7; g[2]=gnm;
	    Ecut[3]=Enm; A[3]=Anm; G[3]=Gnm; a[3]=anm*1.e-7; g[3]=gnm;
	    Ecut[4]=Ene; A[4]=Ane; G[4]=Gne; a[4]=ane*1.e-7; g[4]=gne;
	    Ecut[5]=Ene; A[5]=Ane; G[5]=Gne; a[5]=ane*1.e-7; g[5]=gne;
	    Ecut[6]=Ent; Ecut[7]=Ent;
	    iF = new IntFlux[6];
	    for(i=0; i<6; i++) iF[i] = new IntFlux(type[i], m, f, h0, z0, G[i], prompt, D*1.e-3);
	    iF[0].interpolate(name);
	}

	CteqPDF.tdir=tdir;
	if(!GEN) nT.interpolate(name);

	name+=".gen_";
	if(m==-1) name+="t"; else name+=m+"_"+f+"_"+Output.f(h0);
	if(!PROP && !(Amu==0 && Anm==0 && Ane==0)){
	    name+="_"+Output.f(Ecut[0])+"_"+Output.f(Ecut[2])+"_"+Output.f(Ecut[4]);
	    name+="_"+Output.f(G[0])+"_"+Output.f(G[2])+"_"+Output.f(G[4]);
	    if(prompt!=0) name+="_"+Output.f(prompt);
	}
	name+="_"+Output.f(z0)+"_"+Output.f(D);
	if(!GEN) name+="_"+Output.f(b0)+ename;
	else name+="_go";
	name+="_"+Output.f(mF);
	if(mF>0){
	    if(nlow!=Me) name+="_l"+Output.f(nlow);
	    if(ebig!=bigEnergy) name+="_b"+Output.f(ebig);
	}
	if(sM>0) name+="_s"+sM+"-"+Output.f(Efr);
	if(gE>0) name+="_g"+Output.f(gE)+"-"+Output.f(gD);
	if(Output.raw) name+="_raw"; else name+="_ascii";
	name+=".data";

	boolean flag;
	do {
	    if(Output.texi) return -1;
	    flag=false;
	    try{
		Output.open(name);

		if(!GEN) eM.interpolate();
		if(mF>0){
		    IntFlux.mF=mF;
		    IntFlux.eM=eM;
		    IntFlux.nT=nT;
		}
		if(sM>0){
		    IntFlux.sM=sM;
		    IntFlux.Efr=Efr;
		}
		if(gE>0){
		    IntFlux.gEcut=gE;
		    IntFlux.gD=gD;
		}
		if(!PROP && !(Amu==0 && Anm==0 && Ane==0)){
		    Output.err.print("Parameterizing intflX ... ");
		    for(i=0; i<6; i++) iF[i].interpolate(Ecut[i]);
		    Output.err.println("done");
		    interpolate();
		    Output.err.println("Finished parameterizations");
		}
		Output.close();
	    }catch (mmcException error){
		flag=true;
		Output.delete(name);
	    }
	} while(flag);

	if(!PROP){
	    for(i=0; i<6; i++){
		lept=i; tf[i]=A[i]*getTotFlux(1);
		af[i]=2*Pi*a[i]*Math.pow(Ecut[i<2?i+6:i], 1-g[i])/(g[i]-1);
		tF+=tf[i]+af[i]; sF[i]=tF;
	    }
	    for(i=0; i<6; i++) sF[i]/=tF;
	    tT=nF/(tF*2.e4*Pi*R*R);
	}

	if(f2k){
	    Output.out.println("HI  "+Output.version+"\nHI  "+param);
	    if(!PROP){
		Output.out.println("HI  Total Flux is "+Output.f(tF)+" cm^-2s^-1 ( "+Output.f(tf[0]+tf[1])+" "+Output.f(tf[2]+tf[3])+" "+
				   Output.f(tf[4]+tf[5])+" "+Output.f(af[0]+af[1])+" "+Output.f(af[2]+af[3])+" "+Output.f(af[4]+af[5])+" )");
		Output.out.println("HI  corresponds to the energy cutoffs of Emu="+
				   Output.f(Ecut[0])+" Enm="+Output.f(Ecut[2])+" Ene="+Output.f(Ecut[4])+" Ent="+Output.f(Ecut[6])+" GeV");
		Output.out.println("HI  Amu="+Output.f(A[0])+" Anm="+Output.f(A[2])+" Ane="+Output.f(A[4])
				   +"  Gmu="+Output.f(G[0])+" Gnm="+Output.f(G[2])+" Gne="+Output.f(G[4]));
		Output.out.println("HI  ant="+Output.f(a[0])+" anm="+Output.f(a[2])+" ane="+Output.f(a[4])
				   +"  gnt="+Output.f(g[0])+" gnm="+Output.f(g[2])+" gne="+Output.f(g[4]));
		if(prompt>0) Output.out.println("HI  ratio of prompt to pion muons is "+Output.f(prompt));
		else if(prompt<0) Output.out.println("HI  only prompt muons will be generated: "+Output.f(prompt));
	    }
	    Output.out.println("HI  Detector with radius R="+Output.f(R)+" m is at depth="+Output.f(D)+" m");
	    Output.out.println("HI  Neutrino cross section multiplicative factors are nF="+Output.f(nF)+", mF="+Output.f(mF));
	    if(sM>0) Output.out.println("HI  Solar modulation enabled with parameters: sM="+sM+", Efr="+Output.f(Efr)+" GeV");
	    if(gE>0) Output.out.println("HI  Geomagnetic cutoff enabled with parameters: gE="+Output.f(gE)+
					" GeV, gD="+Output.f(gD)+(sM>0?"":", Efr="+Output.f(Efr)+" GeV"));
	    if(!PROP) Output.out.println("HI  "+Output.f(num)+" showers correspond to the lifetime of "+Output.f(tT*num)+" seconds");
	}

	if(SEED){
	    if(f2k) Output.out.println("HI  random number generator seed is set at "+seed);
	    rand.setSeed(seed);
	}

	return num;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * begins the f2k stream.
     */

    public void beginAtmFlux(){
	if(!GEN) if(f2k){
	    mP.defUSER();
	    tP.defUSER();
	    Output.out.println("USER_DEF almc_e Name Eini Elost Theta");

	    mP.initBufs();
	    tP.initBufs();
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * concludes the f2k stream.
     */

    public void endAtmFlux(){
	if(GEN) Output.err.println("generated "+evt+" events");
	else{
	    Output.err.println("generated "+evt+" events, "+mP.tracks+" muon tracks, "+tP.tracks+" tau tracks");
	    Output.err.println("average "+(mP.vertices+tP.vertices)/evt+" vertices, missed volume: "+(mP.missed+tP.missed));
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * F2k stream lepton generator.
     */

    public static void main(String[] args){
	AtmFlux A = new AtmFlux();
	A.dtct=".atmflux";
	long num=A.initAtmFlux(args, true);

	if(A.PROP){
	    A.mmcf2k();
	}
	else if(num>=0){
	    A.beginAtmFlux();
	    Output.out.println("TBEGIN ? ? ?");
	    for(long i=0; i<num; i++) A.findNext();
	    Output.out.println("TEND ? ? ?");
	    Output.out.println("END");
	    A.endAtmFlux();
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * I3m initializer.
     */

    public void setup(String[] args){
	dtct=".icecube";
	initAtmFlux(args, false);
	I3hist = new Vector(Output.HISTSIZE);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * I3m initializer.
     */

    public void setup(String args){
	setup(Output.splitString(args));
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * I3m function - creates particle set of the next event.
     */

    public Particle[] createNext(){
	if(f2k) return null;
	Particle[] I3p;
	I3hist.clear();
	findNext();
	I3p = new Particle[I3hist.size()];
	for(int i=0; i<I3hist.size(); i++) I3p[i]=(Particle)I3hist.elementAt(i);
	return I3p;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * creates the next event.
     */

    public void findNext(){
	evt++;
	if(!GEN) if(f2k){
	    mP.iniUSER();
	    tP.iniUSER();
	}

	double aux, r, q, fx, fy;
	double ct, st, c1, s1, sa, ca, sp, cp, b, c, tt, tts, ttn;
	double axx, axy, axz, ayx, ayy, ayz, azx, azy, azz, rx, ry, rz;
	double E, x, y, z, t;

	aux=rand.nextDouble();
	lept=0; for(int i=0; i<6; i++) if(aux<sF[i]){ lept=i; break; }

	if(rand.nextDouble()<af[lept]/(af[lept]+tf[lept])){
	    ct=rand.nextDouble();
	    E=Ecut[lept<2?lept+6:lept]*Math.pow(rand.nextDouble(), 1/(1-g[lept]));
	    if(lept<2) lept=lept+6;
	}
	else{
	    ct=getTotFlux(1, rand.nextDouble());
	    E=iF[lept].getE(Ecut[lept], ct, rand.nextDouble());
	}
	th=Math.acos(ct);
	phi=rand.nextDouble()*2*Pi;

	r=R*Math.sqrt(rand.nextDouble());
	q=rand.nextDouble()*2*Pi;
	x=r*Math.cos(q);
	y=r*Math.sin(q);

	st=Math.sqrt(1-ct*ct);  // angle in the detector frame -> angle at the surface
	if(st==0){
	    s1=0;
	    c1=ct;
	}
	else{
	    s1=st*(1-D/R0);
	    c1=Math.sqrt(1-s1*s1);
	}

	x/=c1;

	fx=x*Math.cos(phi)-y*Math.sin(phi);
	fy=x*Math.sin(phi)+y*Math.cos(phi);

	if(rand.nextDouble()>0.5){ ct=-ct; th=Pi-th; }

	sa=st*c1-ct*s1; if(sa>1) sa=1;  // angle difference between the two frames
	ca=ct*c1+st*s1; if(ca>1) ca=1; else if(ca<-1) ca=-1;

	cp=Math.cos(phi);
	sp=Math.sin(phi);

	axx=cp*cp*(ca-1)+1;   // transformation matrix computation
	axy=cp*sp*(ca-1);
	axz=cp*sa;
	ayx=axy;
	ayy=sp*sp*(ca-1)+1;
	ayz=sp*sa;
	azx=-axz;
	azy=-ayz;
	azz=ca;

	rx=R0*sa*cp;          // shift vector computation
	ry=R0*sa*sp;
	rz=R0*(ca-1)+D;       // z is additionally shifted by depth

	x=axx*fx+axy*fy+rx;   // coordinates on the tangent plane in the detector frame
	y=ayx*fx+ayy*fy+ry;
	z=azx*fx+azy*fy+rz;

	cp=Math.cos(phi); sp=Math.sin(phi);

	px=st*cp;             // direction in the CORSIKA frame
	py=st*sp;
	pz=ct;

	b=x*px+y*py+(z+R0-D)*pz;
	c=x*x+y*y+(z+R0-D)*(z+R0-D)-R0*R0;
	r=b-Math.sqrt(b*b-c); // shift along (px,py,pz) to reach the surface from the tangent plane

	x-=px*r;              // final coordinates in the detector frame
	y-=py*r;
	z-=pz*r;
	t=r/(C*1.e-11);

	th*=180/Pi;
	phi*=180/Pi;

	t0+=-Math.log(rand.nextDouble())*tT;
	tt=t0-(x*px+y*py+z*pz)/(C*1.e-2);
	tts=Math.floor(tt*1000)/1000;
	ttn=(tt-tts)*1.e9;
	if(f2k) Output.out.println("EM "+evt+" 1 1970 0 "+Output.f(tts)+" "+Output.f(ttn));

	Name=name[lept];
	eini=E;
	thini=th;

	if(GEN)	putOut(lept, 0, x, y, z, t, -1, E, null);
	else{
	    eM.sett(-px, -py, -pz);

	    if(mF>0 && (lept>=2 && lept<=5)){
		xn=eM.X(eM.ti, eM.tf)*mF;
		boolean nu=(lept==2 || lept==4);
		double tI=(lept==5)?nT.dSdy(E):0;
		tI+=nT.dSdy(E, true, nu);
		tI+=nT.dSdy(E, false, nu);
		xn*=tI;
		if(xn>1) xn=1;
	    }
	    else xn=1;

	    gens=0;
	    elost=0;

	    prop(lept, 0, E, x, y, z, t, eM.ti);
	}

	if(f2k){
	    if(!GEN){
		mP.outUSER();
		tP.outUSER();
		Output.out.println("US almc_e "+Name+" "+Output.f(eini)+" "+
				   Output.f(elost)+" "+Output.f(thini));
	    }
	    Output.out.println("EE");
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * general-purpose lepton propagator.
     */

    public Particle[] propagate(Particle p){
	if(f2k) return null;
	Particle[] I3p;
	I3hist.clear();
	propagate(p.name, p.igen, p.gens, p.x*1.e-2, p.y*1.e-2, p.z*1.e-2, 180-p.theta, p.phi<180?p.phi+180:p.phi-180, p.r*1.e-2, p.e*1.e-3, p.t*1.e9);
	I3p = new Particle[I3hist.size()];
	for(int i=0; i<I3hist.size(); i++) I3p[i]=(Particle)I3hist.elementAt(i);
	return I3p;
    }
    //----------------------------------------------------------------------------------------------------//

    /**
     * Main parser for the f2k file streams.
     */

    public void mmcf2k(){
	Output.err.println(" ---  *** Enter your input in F2000 format now ***  --- ");

	try{
	    LineNumberReader file = new LineNumberReader(new InputStreamReader(Output.in));
	    StringTokenizer st;
	    String line, taux, prnt, type;
	    Vector buffer = new Vector(Output.HISTSIZE);
	    boolean tbegin=false;

	    int i, ipar, igen;
	    long events=0, tracks=0, vertices=0, missed;

	    file.readLine();

	    while((line=file.readLine())!=null){
		st = new StringTokenizer(line);
		if(!st.hasMoreTokens()) continue;
		taux=st.nextToken();
		if("TBEGIN".equals(taux)){
		    beginAtmFlux();
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
		    mP.iniUSER();
		    tP.iniUSER();

		    for(i=0; i<buffer.size(); i++){
			line=(String)buffer.elementAt(i);
			st = new StringTokenizer(line);
			taux=st.nextToken();
			if("TR".equals(taux)){
			    Output.out.println(line);

			    igen=Integer.parseInt(st.nextToken());
			    prnt=st.nextToken();
			    type=st.nextToken();

			    double x, y, z, th, phi, l, e, t;

			    x=Double.parseDouble(st.nextToken());
			    y=Double.parseDouble(st.nextToken());
			    z=Double.parseDouble(st.nextToken());
			    th=Double.parseDouble(st.nextToken());
			    phi=Double.parseDouble(st.nextToken());
			    try{
				l=Double.parseDouble(st.nextToken());
				if(l<0) l=0;
			    }catch(NumberFormatException error){
				l=0;
			    }
			    e=Double.parseDouble(st.nextToken());
			    t=Double.parseDouble(st.nextToken());

			    Name=type;
			    eini=e;
			    thini=th;

			    propagate(type, igen, gens, x, y, z, th, phi, l, e, t);
			    tracks++;

			}
			else if("EM".equals(taux)){
			    events++;
			    Output.out.println(line);
			}
			else if("EE".equals(taux)){
			    mP.outUSER();
			    tP.outUSER();
			    Output.out.println("US almc_e "+Name+" "+Output.f(eini)+" "+
					       Output.f(elost)+" "+Output.f(thini));
			    Output.out.println(line);
			}
			else Output.out.println(line);
		    }

		    gens=0;
		    buffer.clear();
		}

	    }

	    if(events!=0) vertices=(mP.vertices+tP.vertices)/events;
	    missed=mP.missed+tP.missed;
	    Output.err.println("events read "+events+" tracks looped "+tracks+
			       " average vertices "+vertices+" missed volume "+missed+
			       " muon tracks "+mP.tracks+" tau tracks "+tP.tracks);
	}catch(Exception error){
	    Output.err.println("Program finished with exception: "+error.toString());
	    throw new mmcException("input error");
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * general-purpose lepton propagator.
     */

    private void propagate(String type, int igen, int gens, double x, double y, double z, double th, double phi, double r, double E, double t){
	int i;
	double ct, st, cp, sp;
	this.gens=gens;
	this.th=th;
	this.phi=phi;

	for(i=0; i<13; i++) if(name[i].equals(type)) break;
	if(i==13 || i==12 || i==8 || i==9) return;

	ct=Math.cos(th*Pi/180);
	st=Math.sin(th*Pi/180);
	cp=Math.cos(phi*Pi/180);
	sp=Math.sin(phi*Pi/180);

	px=st*cp;
	py=st*sp;
	pz=ct;

	x-=px*r;
	y-=py*r;
	z-=pz*r;
	t+=r/(C*1.e-11);

	eM.sett(-px, -py, -pz);
	elost=0;

	xn=1;
	prop(i, igen, E, x, y, z, t, eM.r(x*1.e-3, y*1.e-3, z*1.e-3));
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * general-purpose lepton propagator.
     */

    private void prop(int lept, int igen, double E, double x, double y, double z, double t, double nf){

	if(nf>eM.tf){
	    putOut(lept, igen, x, y, z, t, 0, E, null);
	    return;
	}

	if(lept>=2 && lept<=7) while(true){
	    if(E<nlow*1.e-3){
		putOut(lept, igen, x, y, z, t, (eM.tf-nf)*1.e3, E, null);
		break;
	    }
	    int mlpt;
	    boolean nu, cc, el;
	    double aux, l, tI, tIc, tIn, tIw=0, rnd, Xt, ml, ni;

	    ni=nf;
	    switch(lept){
	    case 2:  nu=true;  el=false; mlpt=6;  break;
	    case 3:  nu=false; el=false; mlpt=7;  break;
	    case 4:  nu=true;  el=true;  mlpt=4;  break;
	    case 5:  nu=false; el=true;  mlpt=5;  break;
	    case 6:  nu=true;  el=false; mlpt=2;  break;
	    case 7:  nu=false; el=false; mlpt=3;  break;
	    default: nu=false; el=false; mlpt=12;
	    }

	    rnd=rand.nextDouble();
	    tIc=nT.dSdy(E, true, nu, rnd);
	    tIn=nT.dSdy(E, false, nu, rnd);
	    tIw=(lept==5)?nT.dSdy(E, rnd):0;
	    tI=tIc+tIn+tIw;

	    rnd=-xn*Math.log(rand.nextDouble())/(nF*tI);
	    Xt=eM.X(Math.max(ni, eM.ti), eM.tf, rnd);

	    if(rnd>Xt) nf=eM.tf;
	    else{
		nf=eM.t();
		if(Math.abs(nf-eM.tf)<=Math.abs(eM.tf-eM.ti)*computerPrecision) nf=eM.tf;  // computer precision control
	    }
	    l=(nf-ni)*1.e3;

	    putOut(lept, igen, x, y, z, t, l, E, null);
	    igen=gens;

	    if(nf==eM.tf) break;

	    x-=px*l;
	    y-=py*l;
	    z-=pz*l;
	    t+=l/(C*1.e-11);

	    if(OSC) if(!el) if(rand.nextDouble()<nT.Pmt(E, l)) lept=mlpt;

	    switch(lept){
	    case 2:  nu=true;  el=false; mlpt=0;  break;
	    case 3:  nu=false; el=false; mlpt=1;  break;
	    case 4:  nu=true;  el=true;  mlpt=8;  break;
	    case 5:  nu=false; el=true;  mlpt=9;  break;
	    case 6:  nu=true;  el=false; mlpt=10; break;
	    case 7:  nu=false; el=false; mlpt=11; break;
	    default: nu=false; el=false; mlpt=12;
	    }

	    rnd=rand.nextDouble()*tI;
	    if(rnd<tIc){
		aux=nT.e(true, nu);
		E-=aux;
		putOut(12, igen, x, y, z, t, 0, aux, null);
		if(insDet(x, y, z)) elost+=aux;
		prop(mlpt, igen, E, x, y, z, t, nf);
		break;
	    }
	    else if(rnd<tIc+tIn){
		aux=nT.e(false, nu);
		putOut(12, igen, x, y, z, t, 0, aux, null);
		if(insDet(x, y, z)) elost+=aux;
		E-=aux;
	    }
	    else{
		aux=rand.nextDouble();
		if(aux<1/9.){
		    aux=nT.e(Me);
		    putOut(8, igen, x, y, z, t, 0, aux, null);
		    if(insDet(x, y, z)) elost+=aux;
		    E-=aux;
		    lept=5;
		}
		if(aux<2/9.){
		    aux=nT.e(Mmu);
		    prop(0, igen, aux, x, y, z, t, nf);
		    E-=aux;
		    lept=3;
		}
		if(aux<3/9.){
		    aux=nT.e(Mtau);
		    prop(10, igen, aux, x, y, z, t, nf);
		    E-=aux;
		    lept=7;
		}
		else{
		    putOut(12, igen, x, y, z, t, 0, E, null);
		    if(insDet(x, y, z)) elost+=aux;
		    break;
		}
	    }
	}
	else if(lept==0 || lept==1 || lept==10 || lept==11){
	    int nlpt=0;
	    Amanda aP;
	    Particle p[];
	    Particle pi = new Particle(gens+1, gens+1, name[lept], x*1.e2, y*1.e2, z*1.e2, 180-th, phi<180?phi+180:phi-180, E*1.e3, t*1.e-9, 0);
	    aP=lept<2?mP:tP;
	    if(aP.medianum>1) aP.rho[aP.medianum-1]=eM.intRho(nf)/aP.srho[aP.medianum-1];
	    p=aP.propagate(pi); if(f2k) aP.recUSER(); elost+=aP.elost;
	    putOut(lept, igen, x, y, z, t, pi.r*1.e-2, E, pi);
	    igen=gens;
	    for(int i=0; i<p.length; i++){
		boolean flag=true;
		if(p[i].name.equals("nu_e")) nlpt=4;
		else if(p[i].name.equals("~nu_e")) nlpt=5;
		else if(p[i].name.equals("nu_mu")) nlpt=2;
		else if(p[i].name.equals("~nu_mu")) nlpt=3;
		else if(p[i].name.equals("nu_tau")) nlpt=6;
		else if(p[i].name.equals("~nu_tau")) nlpt=7;
		else if(p[i].name.equals("mu") || p[i].name.equals("mu-")) nlpt=0;
		else if(p[i].name.equals("mu+")) nlpt=1;
		else if(p[i].name.equals("tau") || p[i].name.equals("tau-")) nlpt=10;
		else if(p[i].name.equals("tau+")) nlpt=11;
		else flag=false;
		if(flag){
		    double ei=p[i].e*1.e-3, xi=p[i].x*1.e-2, yi=p[i].y*1.e-2, zi=p[i].z*1.e-2;
		    prop(nlpt, igen, ei, xi, yi, zi, p[i].t*1.e9, nf+pi.r*1.e-5);
		    if(lept>=2 && lept<=7) if(insDet(xi, yi, zi)) elost-=ei;
		}
		else lptOut(p[i]);
	    }
	}
	else{
	    putOut(lept, igen, x, y, z, t, 0, E, null);
	    if(insDet(x, y, z)) elost+=E;
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates if the coordinates are inside the detector
     */

    boolean insDet(double x, double y, double z){
	switch(mP.gdet){
	case 1:
	    return Math.abs(x)<mP.length/2 && Math.abs(y)<mP.width/2 && Math.abs(z)<mP.height/2;
	case 2:
	    return x*x+y*y+z*z<mP.radius*mP.radius;
	case 0:
	default:
	    return x*x+y*y<mP.radius*mP.radius && -mP.length/2<z && z<mP.length/2;
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Constructs the lepton from the given information.
     */

    private void putOut(int lept, int igen, double x, double y, double z, double t, double l, double E, Particle p){
	gens++;
	if(f2k) Output.out.println("TR "+gens+(igen>0?" "+igen+" ":" ? ")+name[lept]+" "+Output.f(x)+" "+Output.f(y)+" "+Output.f(z)
				   +" "+Output.f(th)+" "+Output.f(phi)+" "+(l>=0?Output.f(l):"?")+" "+Output.f(E)+" "+Output.f(t));
	else I3hist.add(new Particle(igen, gens, name[lept], x*1.e2, y*1.e2, z*1.e2, 180-th, phi<180?phi+180:phi-180, E*1.e3, t*1.e-9, l*1.e2, p));
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Copies the lepton into the output stream.
     */

    private void lptOut(Particle p){
	gens++;
	if(f2k) Output.out.println("TR "+p.gens+" "+p.igen+" "+p.name+" "+
				   Output.f(p.x*1.e-2)+" "+Output.f(p.y*1.e-2)+" "+Output.f(p.z*1.e-2)+" "+
				   Output.f(180-p.theta)+" "+Output.f(p.phi<180?p.phi+180:p.phi-180)+" "+
				   Output.f(p.r*1.e-2)+" "+Output.f(p.e*1.e-3)+" "+Output.f(p.t*1.e9));
	else I3hist.add(p);
    }

    //----------------------------------------------------------------------------------------------------//

    private double getTotFlux(double x){
	return getTotFlux(x, -1);
    }

    //----------------------------------------------------------------------------------------------------//

    private double getTotFlux(double x, double rnd){
	if(rnd<0){
	    if(jt) return J[lept].interpolate(x);
	    return I.integrateOpened(0, x, this);
	}
	else{
	    if(jt) return Math.min(Math.max(J[lept].findLimit(rnd*J[lept].interpolate(x)), 0), x);
	    I.integrateOpened(0, x, this, rnd);
	    return I.getUpperLimit();
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Function for flux integration over zenith angles - interface to Integral.
     */

    public double function(double x){
	return 2*Pi*iF[lept].getIntFlux(Ecut[lept], x);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Parametrizes the integral of this class.
     */

    public void interpolate(){
	int g=4;  // do not change these settings
	jt=false;
	Output.err.print("Parameterizing totflX ... ");
	J = new Interpolate[6];
	for(int i=0; i<6; i++){
	    lept=i;
	    J[i] = new Interpolate(num3, 0, 1, this, g, true, false, false, g, false, false, true);
	}
	jt=true;
	Output.err.println("done");
    }

    //----------------------------------------------------------------------------------------------------//

    public Interpolate J[];
    public boolean jt=false;

    /**
     * 1d parametrization - interface to Interpolate
     */

    public double functionInt(double x){
	return getTotFlux(x);
    }

}
