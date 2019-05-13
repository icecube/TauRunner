package tfa;
import mmc.*;

/**
 * This class can be used to produce energy losses and spectra of secondaries for all cross sections.
 * Originally thought of to test precision of parameterizations used by mmc (with -<i>test</i>=1 and 2).
 * <ol>
 * <li>
 * -<i>proc</i>=1: calculate energy losses and cross sections (above <i>v</i><sub><sup>cut</sup></sub>)
 * <ul>
 * <li><i>rnd</i><sub><sup>1</sup></sub>&gt;0 sets the random number which selects element of the medium on which the next interaction occurs</li>
 * <li><i>rnd</i><sub><sup>2</sup></sub>&gt;0 sets the random number which selects the energy lost when the next interaction occurs</li>
 * <li><i>-rnd</i><sub><sup>1</sup></sub>&gt;0 sets the initial energy of the muon for the generated table in GeV</li>
 * <li><i>-rnd</i><sub><sup>2</sup></sub>&gt;0 sets the factor by which the energy of the muon is incremented in the generated table</li>
 * </ul>
 * </li>
 * <li>
 * -<i>proc=2</i>: check energy and distance integral calculation
 * <ul>
 * <li><i>rnd</i><sub><sup>1</sup></sub> chooses final energy</li>
 * <li><i>rnd</i><sub><sup>2</sup></sub> chooses final distance</li>
 * </ul>
 * </li>
 * <li>
 * -<i>proc=3</i>: calculate energy spectra for all cross sections (above <i>v</i><sub><sup>cut</sup></sub>)
 * <ul>
 * <li><i>rnd</i><sub><sup>1</sup></sub> sets muon energy in GeV</li>
 * <li><i>rnd</i><sub><sup>2</sup></sub> sets number of points</li>
 * </ul>
 * </ol>
 */

public class Test{

    private static double rnd1=0.5, rnd2=0.5, edef=1.e4;

    public static void main(String[] args){

	double vcut=-1.0;
	double ecut=-1.0;
	String med="Ice";
	String type="mu";
	double elow=PhysicsModel.elow;
	double ebig=PhysicsModel.ebig;

	boolean lpmef=false;
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
	String intr="none";

	int test=1;  // 1 for check_crosS, 2 for check_proP, 3 for plot_crosS

	int bnum=0;
	String param="Test", pbad="";
	boolean pflag;

	for(int n=0; n<args.length; n++){
	    pflag=true;
	    if(args[n].equals("-help") || args[n].equals("-h") || args[n].equals("--help")){
		Output.out.println("\n"+
"This program tests parameterizations used by mmc\n"+
"Available options are: -vcut=[value of vcut: 0<vcut<=1]\n"+
"                       -ecut=[value of ecut in MeV]\n"+
"                       -medi=[medium name]\n"+
"                       -tau  propagate taus instead of muons\n"+
"                       -e    or propagate electrons\n"+
"                       -monopole[=mass in GeV]  or monopoles\n"+
"                       -stau[=mass in GeV]   propagate staus\n"+
"                       -lpm  enable lpm treatment\n"+
"                       -bs=[1-4]  bremsstrahl kkp/abb/ps/csc\n"+
"                       -ph=[1-4]  photonuclear bb/bs/allm/bm\n"+
"                       -bb=[bb/bs:1-4, allm:1-2, bm:1]\n"+
"                       -sh=[1-2]  nuclear structure function\n"+
"                       -c[i/b/p/e/d]=[cross section factor]\n"+
"                       -rho=[multiplicative density factor]\n"+
"                       -elow=[paramet. tables min. in GeV]\n"+
"                       -ebig=[paramet. tables max. in GeV]\n"+
"                       -proc=[TEST to perform: 1, 2, or 3]\n"+
"                       -intr=[interpolate: all, crs or \"\"]\n"+
"                       -romb=[num of interpolation points]\n"+
"                       -RND1=[first  random number]\n"+
"                       -RND2=[second random number]\n"+
"                       -raw  save tables in raw format\n");
		return;
	    }
	    else if(args[n].equals("-tau")){
		type="tau";
	    }
	    else if(args[n].equals("-e")){
		type="e";
	    }
	    else if(args[n].equals("-lpm")){
		lpmef=true;
	    }
	    else if(args[n].equals("-raw")){
		raw=true;
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
	    else if(args[n].startsWith("-proc=")){
		try{
		    test=(int)Double.parseDouble(args[n].substring(6));
		}catch(Exception error){
		    test=0;
		}
	    }
	    else if(args[n].startsWith("-romb=")){
		try{
		    romb=(int)Double.parseDouble(args[n].substring(6));
		}catch(Exception error){
		    romb=5;
		}
	    }
	    else if(args[n].startsWith("-RND1=")){
		try{
		    rnd1=Double.parseDouble(args[n].substring(6));
		}catch(Exception error){
		    rnd1=0.5;
		}
	    }
	    else if(args[n].startsWith("-RND2=")){
		try{
		    rnd2=Double.parseDouble(args[n].substring(6));
		}catch(Exception error){
		    rnd2=0.5;
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
	Propagate p = new Propagate(med , ecut, vcut, type, rho);
	// p.p.low=4.05e2;        // Default settings for the Frejus run are: p.low=105+300 and vcut=0.2
	// p.s.b.lorenzCut=1.e7;  p.s.b.lorenz=true;  // Enable Lorenz invariant violation - cutoff energy in [MeV]
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
	p.interpolate(intr, ".test");  // Comment this out to leave out parameterizations
	Output.DEBUG=false;

	switch(test){

	case 1: check_crosS(p); break;
	case 2: check_proP(p); break;
	case 3: plot_crosS(p); break;
	default: Output.err.println("Error: the proc="+test+" is not defined");

	}
    }

    private static void check_crosS(Propagate p){
	double ind=rnd1;  // random number, chooses element
	double rnd=rnd2;  // .0, ... ,1., chooses lost energy
	double save=0, aux, sum, part, e;
	aux=1/(p.m.Ro*1.e3);
	for(e=ind<0?-ind*1.e3:p.p.low*1.1; e<p.ebig; e*=rnd<0?-rnd:1.2){
	    sum=0;
	    p.p.setEnergy(e);
	    Output.out.println("Energy is "+Output.f(e/1.e3)+" GeV, vCut = "+Output.f(p.m.vCut(p.p.e))+
			       ", medium is "+p.m.name+", energy losses are in [GeV g-1 cm2]");

	    part=p.s.i.c.dEdx()*aux; sum+=part;
	    Output.out.print("Ioniz: "+Output.f(part)+" \t "+Output.f(save=p.s.i.s.dNdx(rnd)));
	    if(save!=0) save=p.s.i.s.e(ind); Output.out.println(" \t "+Output.f(save/e));

	    part=p.s.b.c.dEdx()*aux; sum+=part;
	    Output.out.print("Brems: "+Output.f(part)+" \t "+Output.f(save=p.s.b.s.dNdx(rnd)));
	    if(save!=0) save=p.s.b.s.e(ind); Output.out.println(" \t "+Output.f(save/e));

	    part=p.s.n.c.dEdx()*aux; sum+=part;
	    Output.out.print("Photo: "+Output.f(part)+" \t "+Output.f(save=p.s.n.s.dNdx(rnd)));
	    if(save!=0) save=p.s.n.s.e(ind); Output.out.println(" \t "+Output.f(save/e));

	    part=p.s.e.c.dEdx()*aux; sum+=part;
	    Output.out.print("Epair: "+Output.f(part)+" \t "+Output.f(save=p.s.e.s.dNdx(rnd)));
	    if(save!=0) save=p.s.e.s.e(ind); Output.out.println(" \t "+Output.f(save/e));

	    Output.out.println("Decay: "+Output.f(p.s.d.decay()*p.p.e*aux));

	    Output.out.println("Total: "+Output.f(sum));
	}
    }

    private static void check_proP(Propagate p){
	double rnd=rnd1;  // random number, chooses final energy
	double xfd=rnd2;  // .0, ... ,1., chooses final distance
	double rndMin=0, e, ef=0, xf;
	rnd=-Math.log(rnd);
	for(e=p.p.low*1.1; e<p.ebig; e*=1.2){
	    Output.out.println("Energy is "+Output.f(e/1.e3)+" GeV");
	    Output.out.print(Output.f(rndMin=p.getpr(e, rnd, true))+" \t ");
	    if(rnd<rndMin) Output.out.print(Output.f(ef=p.getef(e, rnd, true))+" \t ");
	    else{ ef=p.p.low; Output.out.print(0+" \t "); }
	    xf=p.s.getdx(e, ef, 0)*xfd;
	    Output.out.print(Output.f(p.s.getdx(e, ef, xf))+" \t ");
	    Output.out.println(Output.f(p.s.getef(e, xf)));
	}
    }

    private static void plot_crosS(Propagate p){
	double emu=rnd1;  // sets muon energy in GeV
	long num=Math.round(rnd2);  // sets number of points

	if(emu<p.p.elow) emu=edef;
	else if(emu>p.ebig) emu=edef;
	if(num<10) num=10;
	else if(num>1000) num=1000;

	int i, n;
	double v, sum, step;
	p.p.setEnergy(emu*1.e3);
	p.s.i.s.setEnergy();

	double alli=p.s.i.s.dNdx();
	double allb=p.s.b.s.dNdx();
	double alle=p.s.e.s.dNdx();
	double alln=p.s.n.s.dNdx();
	double all=alli+allb+alle+alln;

	v=p.m.vCut(p.p.e);
	step=Math.pow(1/v, 1./num);

	Output.out.println("Sum of all cross sections is normalized to 1 from "+Output.f(v*emu)+" to "+Output.f(emu)+" GeV");
	Output.out.println("all = "+Output.f(all)+" alli = "+Output.f(alli)+" allb = "+Output.f(allb)+
			   " alle = "+Output.f(alle)+" alln = "+Output.f(alln));

	for(n=0, v*=Math.sqrt(step); n<num; n++, v*=step){
	    sum=0; if(v>p.s.i.s.vMin && v<p.s.i.s.vMax) sum=p.s.i.s.function(v);
	    Output.out.println("Ioniz: "+Output.f(emu*v)+" \t "+Output.f(sum/all));

	    sum=0; for(i=0; i<p.m.num; i++){ p.s.b.s.setEnergy(i);
	    if(v>p.s.b.s.vMin && v<p.s.b.s.vMax) sum+=p.s.b.s.function(v); }
	    Output.out.println("Brems: "+Output.f(emu*v)+" \t "+Output.f(sum/all));

	    sum=0; for(i=0; i<p.m.num; i++){ p.s.e.s.setEnergy(i);
	    if(v>p.s.e.s.vMin && v<p.s.e.s.vMax) sum+=p.s.e.s.function(v); }
	    Output.out.println("Epair: "+Output.f(emu*v)+" \t "+Output.f(sum/all));

	    sum=0; for(i=0; i<p.m.num; i++){ p.s.n.s.setEnergy(i);
	    if(v>p.s.n.s.vMin && v<p.s.n.s.vMax) sum+=p.s.n.s.function(v); }
	    Output.out.println("Photo: "+Output.f(emu*v)+" \t "+Output.f(sum/all));
	}
    }

}
