package mmc;
import java.util.Random;

/**
 * main mmc class - must be constructed before anything else
 */

public class Propagate extends PhysicsModel{

    public Particle p;
    public Medium m;
    public CrossSections s;
    public Energy2Loss l;
    public Output o;

    public double rho=1;
    public boolean sdec=false;
    public boolean recc=false;
    public boolean exactTime=false;
    public boolean contiCorr=false;
    public boolean molieScat=false;

    public StandardNormal S;
    private Integral I[];
    private boolean pint;
    private double decayS, ionizS, bremsS, epairS, photoS, totalS;

    public boolean dw=false;
    public double rw=0, hw=0;
    public static Random rand = new Random();

    //----------------------------------------------------------------------------------------------------//

    /**
     * initialize all classes necessary for propagation of a muon.
     */

    public Propagate(String w, double ecut, double vcut){
	this(w, ecut, vcut, "mu");
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initialize all classes necessary for propagation of a muon or tau.
     */

    public Propagate(String w, double ecut, double vcut, String type){
	this(w, ecut, vcut, type, 1);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initialize all classes necessary for propagation. String w contains the name of the medium. For the
     * definition of ecut and vcut see help for the class Medium constructor. Set type to "mu" or "tau".
     * The value of rho sets the multiplicative medium density correction factor.
     * To enable parametrization routines, call this.interpolate("all") after this class is created.
     */

    public Propagate(String w, double ecut, double vcut, String type, double rho){
	p = new Particle(this, type);
	m = new Medium(w, ecut, vcut, rho);
	s = new CrossSections(p, m);
	l = new Energy2Loss(s);
	I = new Integral[2];
 	for(int i=0; i<2; i++) I[i] = new Integral(iromb, imaxs, iprec2);
	S = new StandardNormal(iromb, imaxs, iprec);
	o = new Output(p);
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Propagates the particle of initial energy e to the distance r. Returns the final energy if the
     * particle has survived or the track length to the point of disappearance with a minus sign otherwise.
     */

    public double propagateTo(double r, double e){
	int wint;
	boolean DEBUG, flag;
	double ei, ef=0, efd, efi, aux=0, ini, dr;
	double rndd, rndi, rnd1, rnd2, rnd3, rnddMin, rndiMin, rndTot;

	DEBUG=Output.DEBUG;
	if(o.HIST==-1) o.init(p.name);

	ei=e; ef=ei;
	if(r<0) r=0;
	if(e<=p.low || r==0) flag=false;
	else flag=true;
	if(DEBUG) Output.err.println("\nPropagating "+p.name+" of energy "+Output.f(ei)+" MeV to a distance of "+Output.f(r)+" cm");

	while(flag){
	    rndd=-Math.log(rand.nextDouble());
	    rndi=-Math.log(rand.nextDouble());
	    if(DEBUG) Output.err.println("1. solving the tracking integral ...  rndd = "+Output.f(rndd)+"  rndi = "+Output.f(rndi)+" ...  ");
	    rnddMin=p.l<0?0:getpr(ei, rndd, false)/rho;
	    if(DEBUG) Output.err.print(" \t \t \t rnddMin = "+Output.f(rnddMin)+" (d)  ");
	    rndiMin=getpr(ei, rndi, true);
	    if(DEBUG) Output.err.println("rndiMin = "+Output.f(rndiMin)+" (i)");

	    if(DEBUG) Output.err.print("2. evaluating the energy loss ...  ");
	    if(rndd>=rnddMin || rnddMin<=0) efd=p.low;
	    else efd=getef(ei, rndd*rho, false);
	    if(DEBUG) Output.err.print("efd = "+Output.f(efd)+" MeV  ");
	    if(rndi>=rndiMin || rndiMin<=0) efi=p.low;
	    else efi=getef(ei, rndi, true);
	    if(DEBUG) Output.err.println("efi = "+Output.f(efi)+" MeV ...  ");
	    pint=(efi>efd); ef=pint?efi:efd;
	    if(DEBUG) Output.err.println(" \t \t \t lost "+Output.f(ei-ef)+" MeV  ef = "+Output.f(ef)+" MeV");

	    if(DEBUG) Output.err.print("3. calculating the displacement ...  ");
	    dr=s.getdx(ei, ef, rho*(r-p.r))/rho;
	    if(DEBUG) Output.err.println("dr = "+Output.f(dr)+" cm");
	    if(dr<r-p.r){
		if(DEBUG) Output.err.print("4. calculating the local time ...  ");
	    }
	    else{
		dr=r-p.r;
		if(DEBUG) Output.err.print("4. getting the final energy ...  ");
		ef=s.getef(ei, rho*dr);
		if(DEBUG) Output.err.println("lost "+Output.f(ei-ef)+" MeV  ef = "+Output.f(ef)+" MeV");
		if(DEBUG) Output.err.print("5. calculating the local time ...  ");
	    }
	    if(recc) o.output(0, "a"+p.name, ei, dr);
	    p.advance(dr, ei, ef);
	    if(DEBUG) Output.err.println("t = "+Output.f(p.t)+" s");
	    if(Math.abs(r-p.r)<Math.abs(r)*computerPrecision) p.r=r;  // computer precision control
	    if(contiCorr) if(ef!=p.low) ef=S.sndrn(rand.nextDouble(), ef, Math.sqrt(l.e2le.dE2de(ei, ef)), p.low, ei, false);
	    if(recc) o.output(0, "conti", ei-ef, -dr);
	    if(ef==p.low || p.r==r) break;

	    if(DEBUG) Output.err.println("5. choosing the cross section ...");
	    rnd2=rand.nextDouble();
	    rnd3=rand.nextDouble();

	    p.setEnergy(ef);
	    if(pint){
		rnd1=rand.nextDouble();
		if(dw) if(p.r>hw){
		    double exp=Math.abs(rw);
		    double pow=Math.pow(rnd2, exp);
		    rnd2=rw>0?1-pow*rnd2:pow*rnd2;
		    rw=(1+exp)*pow;
		    hw=p.r; dw=false;
		}
		if(DEBUG) decayS=s.d.decay();
		ionizS=s.i.s.dNdx(rnd2);
		bremsS=s.b.s.dNdx(rnd2);
		photoS=s.n.s.dNdx(rnd2);
		epairS=s.e.s.dNdx(rnd2);
		totalS=ionizS+bremsS+photoS+epairS;
		rndTot=rnd1*totalS;
		if(DEBUG) Output.err.println(" . rnd1 = "+Output.f(rnd1)+" rnd2 = "+Output.f(rnd2)+
					     " rnd3 = "+Output.f(rnd3)+" decay = "+Output.f(decayS));
		if(DEBUG) Output.err.println(" . ioniz = "+Output.f(ionizS)+" brems = "+Output.f(bremsS)+
					     " photo = "+Output.f(photoS)+" epair = "+Output.f(epairS));
		if(ionizS>rndTot){ aux=s.i.s.e(rnd3); ef-=aux; wint=2; }
		else if(ionizS+bremsS>rndTot){ aux=s.b.s.e(rnd3); ef-=aux; wint=3; }
		else if(ionizS+bremsS+photoS>rndTot){ aux=s.n.s.e(rnd3); ef-=aux; wint=4; }
		else if(ionizS+bremsS+photoS+epairS>rndTot){ aux=s.e.s.e(rnd3); ef-=aux; wint=5; }
		else{ ei=ef; continue; }  // due to the parameterization of the cross section cutoffs
	    }
	    else{ aux=s.d.e(rnd2, rnd3, p.type==2?rand.nextDouble():0.5, o); ef=0; wint=1; }

	    o.output(wint, wint==1?s.d.out:m.E[s.component], aux, ef);

	    if(ef<=p.low) break;
	    ei=ef;
	}

	if(sdec) if(p.r!=r && ef!=0 && p.l>=0){
	    p.setEnergy(p.m);
	    p.t+=-p.l*Math.log(rand.nextDouble());
	    aux=s.d.e(rand.nextDouble(), 0.5, p.type==2?rand.nextDouble():0.5, o);
	    ef=0;
	    o.output(1, s.d.out, aux, ef);
	}
	p.setEnergy(ef);  // to remember final state of the particle
	o.HIST=-1;  // to make sure user resets particle properties

	if(p.r==r){
	    if(DEBUG) Output.err.println("Particle reached the border with energy ef = "+Output.f(ef)+" MeV");
	    return ef;
	}
	else{
	    if(DEBUG) Output.err.println("Particle "+(p.l<0?"stopped":"disappeared")+" at rf = "+Output.f(p.r)+" cm");
	    return -p.r;
	}
    }

    //----------------------------------------------------------------------------------------------------//

    private double ec, tc;

    /**
     * Propagates the particle of initial energy e to the distance r. Returns the final energy if the
     * particle has survived or the track length to the point of disappearance with a minus sign otherwise.
     * Also calculates particle energy at point rc. Call getPropEc() to get this energy.
     */

    public double propagateTo(double r, double e, double rc){
	int HIST;
	double result;
	HIST=o.HIST;
	ec=propagateTo(rc, e);
	tc=p.t;
	if(ec>0){
	    o.HIST=HIST;
	    result=propagateTo(r, ec);
	}
	else result=ec;
	return result;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * returns the particle energy at rc if the particle has survived or the
     * distance from the point of decay to rc with a minus sign otherwise.
     */

    public double getPropEc(){
	return ec;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * returns the particle time at rc if the particle has survived or at the point of decay otherwise.
     */

    public double getPropTc(){
	return tc;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * function for energy range calculation - interface to Integral
     */

    public double function(double E){
	final boolean DEBUG=false;
	double aux;
	aux=s.function(E);
	if(!pint){
	    decayS=s.d.decay();
	    if(DEBUG) Output.err.print(" + "+Output.f(p.e));
	    return aux*decayS;
	}
	else{
	    ionizS=s.i.s.dNdx();
	    if(DEBUG) Output.err.print(" \t "+Output.f(ionizS));
	    bremsS=s.b.s.dNdx();
	    if(DEBUG) Output.err.print(" \t "+Output.f(bremsS));
	    photoS=s.n.s.dNdx();
	    if(DEBUG) Output.err.print(" \t "+Output.f(photoS));
	    epairS=s.e.s.dNdx();
	    if(DEBUG) Output.err.println(" \t "+Output.f(epairS));
	    totalS=ionizS+bremsS+photoS+epairS;
	    return aux*totalS;
	}
    }

    //----------------------------------------------------------------------------------------------------//

    public static int g=5;

    /**
     * call this routine to enable interpolations. To enable everything, set w="all"
     */

    public void interpolate(String w){
	int i;
	double e_hi=ebig;
	double e_low=p.low;

	Output.err.println("Parameterizations apply in the energy range from "+Output.f(e_low)+" MeV to "+Output.f(e_hi)+" MeV");
	if(w.indexOf("all")!=-1){
	    w+=" crs trackE trackX";
	    if(exactTime) w+=" trackT";
	    if(contiCorr) w+=" contiR";
	    if(molieScat) w+=" molieS gaussR";
	}
	if(w.indexOf("crs")!=-1){
	    w+=" ionizE ionizS bremsE bremsS photoE photoS epairE epairP epairS";
	    if(s.n.form>2) w+=" photoP";
	}
	if(w.indexOf("contiR")!=-1) w+=" gaussR en2ldX en2ldE";
	w=w.toLowerCase();

	s.i.c.jt=false;
	if(w.indexOf("ionize")!=-1 && s.ci>0){
	    Output.err.print("Parameterizing ionizE ... ");
	    s.i.c.J = new Interpolate(num1, e_low, e_hi, s.i.c, g, true, false, true, g, false, false, true);
	    s.i.c.jt=true;
	    Output.err.println("done");
	}

	s.i.s.jt=false;
	if(w.indexOf("ionizs")!=-1 && s.ci>0){
	    Output.err.print("Parameterizing ionizS ... ");
	    s.i.s.J = new Interpolate(num1, e_low, e_hi, num1, 0, 1, s.i.s, g, false, false, true, g, false, false, false, g, true, false, false);
	    s.i.s.Jo = new Interpolate(num1, e_low, e_hi, s.i.s, g, false, false, true, g, true, false, false);
	    s.i.s.jt=true;
	    Output.err.println("done");
	}

	s.b.c.jt=false;
	if(w.indexOf("bremse")!=-1 && s.cb>0){
	    Output.err.print("Parameterizing bremsE ... ");
	    s.b.c.J = new Interpolate(num1, e_low, e_hi, s.b.c, g, true, false, true, g, false, false, false);
	    s.b.c.jt=true;
	    Output.err.println("done");
	}

	s.b.s.jt=false;
	if(w.indexOf("bremss")!=-1 && s.cb>0){
	    Output.err.print("Parameterizing bremsS ... ");
	    s.b.s.J = new Interpolate[m.num];
	    s.b.s.Jo = new Interpolate[m.num];
	    for(i=0; i<m.num; i++){
		s.component=i;
		s.b.s.J[i] = new Interpolate(num1, e_low, e_hi, num1, 0, 1, s.b.s, g, false, false, true, g, false, false, false, g, true, false, false);
		s.b.s.Jo[i] = new Interpolate(num1, e_low, e_hi, s.b.s, g, false, false, true, g, true, false, false);
	    }
	    s.b.s.jt=true;
	    Output.err.println("done");
	}

	s.n.jt=false;
	if(w.indexOf("photop")!=-1 && s.cp>0){
	    Output.err.print("Parameterizing photoP ... ");
	    s.n.J = new Interpolate[m.num];
	    for(i=0; i<m.num; i++){
		s.component=i;
		s.n.J[i] = new Interpolate(num1, e_low, e_hi, num1, 0., 1., s.n, g, false, false, true, g, false, false, false, g, false, false, false);
	    }
	    s.n.jt=true;
	    Output.err.println("done");
	}

	s.n.c.jt=false;
	if(w.indexOf("photoe")!=-1 && s.cp>0){
	    Output.err.print("Parameterizing photoE ... ");
	    s.n.c.J = new Interpolate(num1, e_low, e_hi, s.n.c, g, true, false, true, g, false, false, false);
	    s.n.c.jt=true;
	    Output.err.println("done");
	}

	s.n.s.jt=false;
	if(w.indexOf("photos")!=-1 && s.cp>0){
	    Output.err.print("Parameterizing photoS ... ");
	    s.n.s.J = new Interpolate[m.num];
	    s.n.s.Jo = new Interpolate[m.num];
	    for(i=0; i<m.num; i++){
		s.component=i;
		s.n.s.J[i] = new Interpolate(num1, e_low, e_hi, num1, 0, 1, s.n.s, g, false, false, true, g, false, false, false, g, true, false, false);
		s.n.s.Jo[i] = new Interpolate(num1, e_low, e_hi, s.n.s, g, false, false, true, g, true, false, false);
	    }
	    s.n.s.jt=true;
	    Output.err.println("done");
	}

	s.e.jt=false;
	if(w.indexOf("epairp")!=-1 && s.ce>0){
	    Output.err.print("Parameterizing epairP ... ");
	    s.e.J = new Interpolate[m.num];
	    for(i=0; i<m.num; i++){
		s.component=i;
		s.e.J[i] = new Interpolate(num1, e_low, e_hi, num1, 0., 1., s.e, g, false, false, true, g, false, false, false, g, false, false, false);
	    }
	    s.e.jt=true;
	    Output.err.println("done");
	}

	s.e.c.jt=false;
	if(w.indexOf("epaire")!=-1 && s.ce>0){
	    Output.err.print("Parameterizing epairE ... ");
	    s.e.c.J = new Interpolate(num1, e_low, e_hi, s.e.c, g, true, false, true, g, false, false, false);
	    s.e.c.jt=true;
	    Output.err.println("done");
	}

	s.e.s.jt=false;
	if(w.indexOf("epairs")!=-1 && s.ce>0){
	    Output.err.print("Parameterizing epairS ... ");
	    s.e.s.J = new Interpolate[m.num];
	    s.e.s.Jo = new Interpolate[m.num];
	    for(i=0; i<m.num; i++){
		s.component=i;
		s.e.s.J[i] = new Interpolate(num1, e_low, e_hi, num1, 0, 1, s.e.s, g, false, false, true, g, false, false, false, g, true, false, false);
		s.e.s.Jo[i] = new Interpolate(num1, e_low, e_hi, s.e.s, g, false, false, true, g, true, false, false);
	    }
	    s.e.s.jt=true;
	    Output.err.println("done");
	}

	this.jt=false;
	if(w.indexOf("tracke")!=-1){
	    Output.err.print("Parameterizing trackE ... ");
	    this.J = new Interpolate[2];
	    this.Jdf = new Interpolate[2];
	    pint=true; if(Math.abs(-I[1].integrateWithLog(p.low, p.low*10, this))<Math.abs(-I[1].integrateWithLog(ebig, ebig/10, this))) up=true; else up=false;
	    for(i=0; i<2; i++){
		pint=(i==1);
		this.df=false; this.J[i] = new Interpolate(num3, e_low, e_hi, this, g, false, false, true, g, false, false, false);
		this.df=true; this.Jdf[i] = new Interpolate(num3, e_low, e_hi, this, g, false, false, true, g, false, false, false);
	    }
	    this.jt=true;
	    for(i=0; i<2; i++) if(!(up&&(i==1))) bigLow[i]=J[i].interpolate(p.low);
	    if(up) Output.err.println("done");
	    else Output.err.println("down");
	}

	s.jt=false;
	if(w.indexOf("trackx")!=-1){
	    Output.err.print("Parameterizing trackX ... ");
	    s.df=false; s.J = new Interpolate(num3, e_low, e_hi, s, g, false, false, true, g, false, false, false);
	    s.df=true; s.Jdf = new Interpolate(num3, e_low, e_hi, s, g, false, false, true, g, false, false, false);
	    s.jt=true;
	    Output.err.println("done");
	}

	p.jt=false;
	if(w.indexOf("trackt")!=-1){
	    Output.err.print("Parameterizing trackT ... ");
	    p.df=false; p.J = new Interpolate(num3, e_low, e_hi, p, g, false, false, true, g, false, false, false);
	    p.df=true; p.Jdf = new Interpolate(num3, e_low, e_hi, p, g, false, false, true, g, false, false, false);
	    p.jt=true;
	    Output.err.println("done");
	}

	p.s.jt=false;
	if(w.indexOf("molies")!=-1){
	    Output.err.print("Parameterizing molieS ... ");
	    p.s.df=false; p.s.J = new Interpolate(num2, e_low, e_hi, p.s, g, false, false, true, g, false, false, false);
	    p.s.df=true; p.s.Jdf = new Interpolate(num2, e_low, e_hi, p.s, g, false, false, true, g, false, false, false);
	    p.s.jt=true;
	    Output.err.println("done");
	}

	S.jt=false;
	if(w.indexOf("gaussr")!=-1){
	    Output.err.print("Parameterizing gaussR ... ");
	    S.J = new Interpolate(num2, -5, 5, S, g, true, false, false, g, true, false, false);
	    S.jt=true;
	    Output.err.println("done");
	}

	l.e2lx.jt=false;
	if(w.indexOf("en2ldx")!=-1){
	    Output.err.print("Parameterizing en2ldX ... ");
	    l.e2lx.J = new Interpolate(num2, e_low, e_hi, l.e2lx, g, false, false, true, g, false, false, false);
	    l.e2lx.jt=true;
	    Output.err.println("done");
	}

	l.e2le.jt=false;
	if(w.indexOf("en2lde")!=-1){
	    Output.err.print("Parameterizing en2ldE ... ");
	    l.e2le.df=false; l.e2le.J = new Interpolate(num2, e_low, e_hi, l.e2le, g, false, false, true, g, false, false, false);
	    l.e2le.df=true; l.e2le.Jdf = new Interpolate(num2, e_low, e_hi, l.e2le, g, false, false, true, g, false, false, false);
	    l.e2le.jt=true;
	    Output.err.println("done");
	}

	Output.err.println("Finished parameterizations");
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * call this routine to enable interpolations and their save and reread. To enable everything, set w="all"
     */

    public void interpolate(String w, String filename){
	String name;
	boolean flag;
	name=(p.name+"_"+(m.name+"_"+w).toLowerCase()).replace(' ', '-');
	name=filename+"."+name+"_"+Output.f(m.ecut)+"_"+Output.f(m.vcut)+"_";
	name+=s.b.form;
	name+=s.n.form;
	name+=s.n.bb;
	name+=s.n.shadow;
	if(s.lpm) name+="t"; else name+="f";
	if(exactTime) name+="t"; else name+="f";
	if(contiCorr) name+="t"; else name+="f";
	if(molieScat) name+="t"; else name+="f";
	if(elow==p.low) name+="_l"+Output.f(elow);
	if(ebig!=bigEnergy) name+="_b"+Output.f(ebig);
	if(s.ci!=1 || s.cb!=1 || s.cp!=1 || s.ce!=1 || s.cd!=1 || m.rho!=1){
	    name+="_"+Output.f(s.ci);
	    name+=","+Output.f(s.cb);
	    name+=","+Output.f(s.cp);
	    name+=","+Output.f(s.ce);
	    name+=","+Output.f(s.cd);
	    name+=","+Output.f(m.rho);
	}
	if(Output.raw) name+="_raw"; else name+="_ascii";
	name+=".data";
	do {
	    if(Output.texi) return;
	    flag=false;
	    try{
		Output.open(name);
		interpolate(w);
		Output.close();
	    }catch (mmcException error){
		flag=true;
		Output.delete(name);
	    }
	} while(flag);
    }

    //----------------------------------------------------------------------------------------------------//

    private boolean up=true;
    private double bigLow[] = new double[2], storeDif[] = new double[2];

    /**
     * returns the value of the tracking integral
     */

    public double getpr(double ei, double rnd, boolean pint){
	if(jt){
	    storeDif[pint?1:0]=J[pint?1:0].interpolate(ei);
	    if(up&&pint) return Math.max(storeDif[pint?1:0], 0);
	    else return Math.max(bigLow[pint?1:0]-storeDif[pint?1:0], 0);
	}
	else{
	    this.pint=pint;
	    return I[pint?1:0].integrateWithLog(ei, p.low, this, -rnd);
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * final energy, corresponding to the given value rnd of the tracking integral
     */

    public double getef(double ei, double rnd, boolean pint){
	if(jt){
	    if(Math.abs(rnd)>Math.abs(storeDif[pint?1:0])*halfPrecision){
		double aux;
		if(up&&pint) aux=J[pint?1:0].findLimit(storeDif[pint?1:0]-rnd);
		else aux=J[pint?1:0].findLimit(storeDif[pint?1:0]+rnd);
		if(Math.abs(ei-aux)>Math.abs(ei)*halfPrecision) return Math.min(Math.max(aux, p.low), ei);
	    }
	    return Math.min(Math.max(ei+rnd/Jdf[pint?1:0].interpolate(ei+rnd/(2*Jdf[pint?1:0].interpolate(ei))), p.low), ei);
	}
	else{
	    this.pint=pint;
	    return I[pint?1:0].getUpperLimit();
	}
    }

    //----------------------------------------------------------------------------------------------------//

    boolean df=false;
    public Interpolate J[];
    public Interpolate Jdf[];
    public boolean jt=false;

    /**
     * 1d parametrization - interface to Interpolate
     */

    public double functionInt(double e){
	if(df) return function(e);
	else{
	    if(up&&pint) return I[pint?1:0].integrateWithLog(e, p.low, this);
	    else return -I[pint?1:0].integrateWithLog(e, ebig, this);
	}
    }

}
