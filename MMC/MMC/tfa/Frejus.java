package tfa;
import mmc.*;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.util.StringTokenizer;

/**
 * This class evaluates the energy of muons/taus after passing through a given distance or their range.
 * Reads from and outputs to standard input/output. Originally thought of as a front-end for Frejus experiment.
 */

public class Frejus{

    public static void main(String[] args){

	double vcut=-1.0;
	double ecut=-1.0;
	String med="Frejus Rock";
	String type="mu";
	double elow=PhysicsModel.elow;
	double ebig=PhysicsModel.ebig;
	boolean SDEC=false;
	boolean RECC=false;
	boolean timef=false;
	boolean conti=false;
	boolean lpmef=false;
	boolean debug=false;

	int bspar=1;
	int pncrs=1;
	int pncbb=1;
	int pncsh=1;
	double crsci=1;
	double crscb=1;
	double crscp=1;
	double crsce=1;
	double crscd=1;
	double rho=1;
	int romb=5;
	boolean raw=false;
	String intr="all";

	int bnum=0;
	String param="Frejus", pbad="";
	boolean pflag;

	for(int n=0; n<args.length; n++){
	    pflag=true;
	    if(args[n].equals("-help") || args[n].equals("-h") || args[n].equals("--help")){
		Output.out.println("\n"+
"This program propagates muons of given energy through a given distance\n"+
"Available options are: -vcut=[value of vcut: 0<vcut<=1]\n"+
"                       -ecut=[value of ecut in MeV]\n"+
"                       -medi=[medium name]\n"+
"                       -tau  propagate taus instead of muons\n"+
"                       -e    propagate electrons instead of muons\n"+
"                       -monopole[=mass in GeV] propagate monopoles\n"+
"                       -stau[=mass in GeV]     propagate staus\n"+
"                       -sdec enable stopped muon decay treatment\n"+
"                       -recc enable printout of continuous energy losses\n"+
"                       -time also report time of flight\n"+
"                       -cont enable continuous loss randomization\n"+
"                       -lpm  enable lpm treatment\n"+
"                       -bs=[1-4]  bremsstrahlung: kkp, abb, ps, csc\n"+
"                       -ph=[1-4]  photonuclear: bb, bb+bs, allm, bm\n"+
"                       -bb=[bb/bs:1-4 3|4:bb|zeus, allm:1-2(91/7), bm:1]\n"+
"                       -sh=[1-2]  nuclear structure function: dutt/butk\n"+
"                       -c[i/b/p/e/d]=[cross section modifier, 0:disable]\n"+
"                       -rho=[multiplicative factor] for medium density\n"+
"                       -elow=[muon energy in GeV below which it is lost]\n"+
"                       -ebig=[upper bound in GeV of the paramet. tables]\n"+
"                       -intr=[interpolate: \"all\", \"crs\" or \"\"]\n"+
"                       -romb=[number of interpolation points]\n"+
"                       -raw  save tables in raw format\n"+
"                       -debug print debugging information to stderr\n");
		return;
	    }
	    else if(args[n].equals("-tau")){
		type="tau";
	    }
	    else if(args[n].equals("-e")){
		type="e";
	    }
	    else if(args[n].equals("-sdec")){
		SDEC=true;
	    }
	    else if(args[n].equals("-recc")){
		RECC=true;
	    }
	    else if(args[n].equals("-time")){
		timef=true;
	    }
	    else if(args[n].equals("-cont")){
		conti=true;
	    }
	    else if(args[n].equals("-lpm")){
		lpmef=true;
	    }
	    else if(args[n].equals("-raw")){
		raw=true;
	    }
	    else if(args[n].equals("-debug")){
		debug=true;
	    }
	    else if(args[n].startsWith("-monopole")){
		type=args[n].substring(1);
	    }
	    else if(args[n].startsWith("-stau")){
		type=args[n].substring(1);
	    }
	    else if(args[n].startsWith("-vcut=")){
		try{
		    vcut=Double.parseDouble(args[n].substring(6));
		}catch(Exception error){
		    vcut=-1;
		}
	    }
	    else if(args[n].startsWith("-ecut=")){
		try{
		    ecut=Double.parseDouble(args[n].substring(6));
		}catch(Exception error){
		    ecut=-1;
		}
	    }
	    else if(args[n].startsWith("-bs=")){
		try{
		    bspar=(int)Double.parseDouble(args[n].substring(4));
		}catch(Exception error){
		    bspar=1;
		}
	    }
	    else if(args[n].startsWith("-ph=")){
		try{
		    pncrs=(int)Double.parseDouble(args[n].substring(4));
		}catch(Exception error){
		    pncrs=1;
		}
	    }
	    else if(args[n].startsWith("-bb=")){
		try{
		    pncbb=(int)Double.parseDouble(args[n].substring(4));
		}catch(Exception error){
		    pncbb=1;
		}
	    }
	    else if(args[n].startsWith("-sh=")){
		try{
		    pncsh=(int)Double.parseDouble(args[n].substring(4));
		}catch(Exception error){
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
		    rho=Double.parseDouble(args[n].substring(5));
		}catch(NumberFormatException error){
		    rho=1;
		}
	    }
	    else if(args[n].startsWith("-elow=")){
		try{
		    elow=Double.parseDouble(args[n].substring(6))*1.e3;
		}catch(Exception error){
		    elow=PhysicsModel.elow;
		}
	    }
	    else if(args[n].startsWith("-ebig=")){
		try{
		    ebig=Double.parseDouble(args[n].substring(6))*1.e3;
		}catch(Exception error){
		    ebig=PhysicsModel.ebig;
		}
	    }
	    else if(args[n].startsWith("-romb=")){
		try{
		    romb=(int)Double.parseDouble(args[n].substring(6));
		}catch(Exception error){
		    romb=5;
		}
	    }
	    else if(args[n].startsWith("-medi=")){
		med=args[n].substring(6).replace('-', ' ');
	    }
	    else if(args[n].startsWith("-intr=")){
		intr=args[n].substring(6).replace('-', ' ');
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
	    pbad=bnum==1?pbad+" is":"s"+pbad+" are";
	    Output.err.println("Warning: Parameter"+pbad+" not recognized");
	}

	if(bspar<1 || bspar>4){
	    Output.err.println("Warning: bs is not a valid number");
	    bspar=1;
	}
	if(pncrs<1 || pncrs>4 || (pncrs==2 && !(type.equals("mu") || type.equals("tau")))){
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

	PhysicsModel.elow=elow;
	PhysicsModel.ebig=ebig;
	Propagate p = new Propagate(med, ecut, vcut, type, rho);
	Output.DEBUG=debug;
	p.sdec=SDEC;           // To enable stopped muon decay
	p.recc=RECC;           // To enable printout of continuous energy losses
	p.contiCorr=conti;     // To randomize the continuous energy losses, use only for small vcut
	p.exactTime=timef;     // To compute local time of the particle exactly
	p.s.lpm=lpmef;         // Enable lpm and dielectric suppression effects
	p.s.b.form=bspar;      // Choose parametrization of the bremsstrahlung cross section
	p.s.n.form=pncrs;      // Choose parametrization of the photon-nucleon cross section
	p.s.n.bb=pncbb;        // Choose parametrization of the photon-nucleon cross section
	p.s.n.shadow=pncsh;    // Choose parametrization of the photon-nucleon cross section
	p.s.ci=crsci;          // Cross section multiplicative modifier
	p.s.cb=crscb;          // Cross section multiplicative modifier
	p.s.cp=crscp;          // Cross section multiplicative modifier
	p.s.ce=crsce;          // Cross section multiplicative modifier
	p.s.cd=crscd;          // Cross section multiplicative modifier
	Propagate.g=romb;
	Output.raw=raw;
	p.interpolate(intr, ".frejus");
	Output.err.println("Enter the following: energy in [GeV] and distance to travel x in [m].");
	Output.err.println("The last entry in the output is the final energy in [GeV] (if > 0)");
	Output.err.println("or distance traveled to the point of disappearance in [m] (if < 0).");
	Output.err.println("          ---  *** Enter your parameters now ***  ---");

	try{

	    LineNumberReader file = new LineNumberReader(new InputStreamReader(Output.in));
	    String line;
	    StringTokenizer t;
	    double e, x, result;

	    while((line=file.readLine())!=null){
		t = new StringTokenizer(line);
		if(!t.hasMoreTokens()) continue;
		e=Double.parseDouble(t.nextToken());
		x=Double.parseDouble(t.nextToken());
		result=p.propagateTo(x*1.e2, e*1.e3);
		line+=" "+Output.f(result>0?result*1.e-3:result*1.e-2);
		if(timef) line+=" "+Output.f(p.p.t);
		Output.out.println(line);
	    }

	}catch(Exception error){
	    Output.err.println("Program finished with exception: "+error.toString());
	    throw new mmcException("input error");
	}
    }

}
