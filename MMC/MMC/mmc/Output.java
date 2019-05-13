package mmc;
import java.util.Vector;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.ByteArrayOutputStream;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.io.File;
import java.io.Reader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.URL;

/**
 * Class contains all parametrization table io routines, history/debugging info/secondary
 * output definitions, and a fast double precision variable formatting routine
 */

public class Output{

    final public static String version="Muon Propagation Code in Java v. 1.6.0";

    final public static int HISTSIZE=1000;

    public static InputStream in = System.in;
    public static PrintStream out = System.out;
    public static PrintStream err = System.err;

    public static boolean DEBUG=false;
    public static boolean AMASIM=false;
    public StringBuffer history = new StringBuffer(HISTSIZE);
    public int igen, gens;

    public static boolean I3flag=false;
    public static boolean RecDec=false;
    public Vector I3hist = new Vector(HISTSIZE);

    protected int HIST=-1;
    protected Particle p;

    //----------------------------------------------------------------------------------------------------//

    /**
     * initialize all classes necessary for propagation
     */

    public Output(Particle p){
	this.p = p;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * enables redirection of the stderr output into a file
     */

    public static boolean setStderr(String name){
	try{
	    if(err == System.err) err = new PrintStream(new FileOutputStream(name));
	    return true;
	}
	catch(Exception e){
	    Output.err.println("Cannot redirect stderr to the file "+name);
	    return false;
	}
    }

    //----------------------------------------------------------------------------------------------------//

    ByteArrayOutputStream outstr;

    /**
     * enables redirection of the stdout output into a string
     */

    public void setStdout(){
	outstr = new ByteArrayOutputStream();
	out = new PrintStream(outstr);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * gets the stdout output string, resets the stdout output
     */

    public String getStdout(){
	out = System.out;	
	return outstr.toString();
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initialize for no history output - no need to call this function
     */

    public void init(String type){
	HIST=0;
	p.location(type, 0, 0, 0, 0, 0, 0);
	history.setLength(0);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initialize for regular history output
     */

    public void initDefault(String type){
	HIST=1;
	p.location(type, 0, 0, 0, 0, 0, 0);
	history.setLength(0);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initialize for regular history output
     */

    public void initDefault(int igen, int gens, String type, double time, double x, double y, double z, double theta, double phi){
	HIST=3;
	this.igen=igen;
	this.gens=gens;
	p.location(type, time, x, y, z, theta, phi);
	history.setLength(0);
	if(RecDec) if(I3flag) I3hist.clear();
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initialize for F2000 history output
     */

    public void initF2000(int igen, int gens, String type, double time, double x, double y, double z, double theta, double phi){
	HIST=2;
	this.igen=igen;
	this.gens=gens;
	p.location(type, time, x, y, z, theta, phi);
	history.setLength(0);
	if(I3flag) I3hist.clear();
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * print output information in the format specified by HIST
     */

    protected void output(int wint, String comp, double de, double ef){
	String str;
	if(DEBUG || (HIST==1)){
	    str=" . ";
	    if(wint==0){
		if(comp=="conti") str+="lost continuously "+f(de)+" MeV over "+f(ef)+" cm up to "+f(p.r)+" cm\n";
		else str+="muon with energy "+f(de)+" MeV traveled  "+f(ef)+" cm\n";
	    }
	    else{
		if(wint==1) str+="decayed into "+comp;
		else if(wint==2) str+="ionized by "+comp;
		else if(wint==3) str+="bremsed by "+comp;
		else if(wint==4) str+="photoed by "+comp;
		else if(wint==5) str+="epaired by "+comp;
		str+=" ... lost "+f(de)+" MeV, ef = "+f(ef)+" MeV, xf = "+f(p.r)+" cm\n";
	    }
	    if(HIST==1) history.append(str);
	    if(DEBUG) Output.err.print(str);
	}
	if(HIST==2 || (RecDec && wint==1)){
	    if(wint==0 || wint==1) str=comp;
	    else if(wint==2) str="delta";
	    else if(wint==3) str="brems";
	    else if(wint==4) str="munu";
	    else if(wint==5) str="epair";
	    else str="";
	    if(str!=""){
		gens++;
		if(I3flag) I3hist.addElement(new Particle(igen, gens, str, p.x, p.y, p.z, p.theta, p.phi, de, p.t, wint==0?ef:0));
		else{
		    str="TR "+gens+" "+igen+" "+str+" "+f(p.x*1.e-2)+" "+f(p.y*1.e-2)+" "+f(p.z*1.e-2);
		    str+=" "+f(180-p.theta)+" "+f(p.phi<180?p.phi+180:p.phi-180)+" "+f(wint==0?ef*1.e-2:0)+" "+f(de*1.e-3)+" "+f(p.t*1.e9)+"\n";
		    history.append(str);
		}
	    }
	}
    }

    //----------------------------------------------------------------------------------------------------//

    protected static boolean inf=false;
    protected static boolean outf=false;

    protected static LineNumberReader ina=null;
    protected static PrintWriter outa=null;

    protected static ObjectInputStream inr=null;
    protected static ObjectOutputStream outr=null;

    public static boolean raw=false;
    public static boolean web=false;
    public static boolean sta=false;
    public static boolean texi=false;

    //----------------------------------------------------------------------------------------------------//

    /**
     * read one number from the table
     */

    protected static double read(){
	try{
	    if(raw) return inr.readDouble();
	    else{
		String str=ina.readLine();
		if(str==null) throw new mmcException("parametrization table premature eof");
		return Double.parseDouble(str);
	    }
	}catch(Exception error){
	    Output.err.println("***WARNING***: read failed: "+error.toString());
	    close();
	    throw new mmcException("parametrization table read error");
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * write one number to the table
     */

    protected static void write(double g){
	try{
	    if(raw) outr.writeDouble(g);
	    else outa.println(g);
	}catch(Exception error){
	    Output.err.println("***WARNING***: write failed: "+error.toString());
	    close();
	}
    }

    //----------------------------------------------------------------------------------------------------//

    protected static FileInputStream fileIr=null;
    protected static FileReader fileIa=null;
    protected static FileOutputStream fileOr=null;
    protected static FileWriter fileOa=null;

    /**
     * See if the file "name" exists
     */

    public static boolean exists(String name){
	sta=name.startsWith("http://");
	if(web || sta){
	    try{
		(new URL(name)).openStream();
		return true;
	    } catch (Exception error) {
		return false;
	    }
	}
	return (new File(name)).exists();
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * open the table
     */

    public static void open(String name){
	sta=name.startsWith("http://");
	inf=false; outf=false;
	try{
	    if(web || sta){
		if(name.equals(badWebFile)) throw new mmcException("web file overwrite failed");
		if(exists(name)){
		    if(raw) inr = new ObjectInputStream((new URL(name)).openStream());
		    else throw new mmcException("cannot use this format for web transfers");
		    inf=true;
		    Output.err.println("Parametrization tables will be read in from the file "+name);
		}
		else throw new mmcException("cannot save file at this location");
	    }
	    else{
		if(exists(name)){
		    if(raw) inr = new ObjectInputStream(fileIr = new FileInputStream(name));
		    else ina = new LineNumberReader(fileIa = new FileReader(name));
		    inf=true;
		    Output.err.println("Parametrization tables will be read in from the file "+name);
		}
		else{
		    if(raw) outr = new ObjectOutputStream(fileOr = new FileOutputStream(name));
		    else outa = new PrintWriter(fileOa = new FileWriter(name));
		    outf=true;
		    Output.err.println("Parametrization tables will be saved to the file "+name);
		}
	    }
	}catch(Exception error){
	    Output.err.println("File "+name+" cannot be used: "+error.toString());
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * close the table
     */

    public static void close(){
	inf=false; outf=false;
	try{
	    if(web || sta){ }
	    else{
		if(raw){
		    if(inr!=null) inr.close();
		    if(fileIr!=null) fileIr.close();
		    if(outr!=null) outr.close();
		    if(fileOr!=null) fileOr.close();
		}
		else{
		    if(ina!=null) ina.close();
		    if(fileIa!=null) fileIa.close();
		    if(outa!=null) outa.close();
		    if(fileOa!=null) fileOa.close();
		}
	    }
	}catch(Exception error){
	    Output.err.println("Failed to close file: "+error.toString());
	}
	inr=null; outr=null;
	ina=null; outa=null;
	fileIr=null; fileOr=null;
	fileIa=null; fileOa=null;
    }

    //----------------------------------------------------------------------------------------------------//

    private static String badWebFile="";

    /**
     * delete the file
     */

    public static void delete(String name){
	sta=name.startsWith("http://");
	inf=false; outf=false;
	try{
	    Output.err.println("... Deleting corrupt file "+name);
	    if(web || sta) badWebFile=name;
	    (new File(name)).delete();
	}catch(Exception error){
	    Output.err.println("Failed to delete file "+name+": "+error.toString());
	    return;
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * open Reader stream
     */

    public static Reader reader(String name) throws Exception{
	sta=name.startsWith("http://");
	return web || sta?new InputStreamReader((new URL(name)).openConnection().getInputStream()):new FileReader(name);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Splits the string into the array
     */

    public static String[] splitString(String args){
	int i=0, j;
	Vector parm = new Vector(Output.HISTSIZE);
	while((j=args.indexOf(" ", i))!=-1){
	    if(j>0) parm.addElement(args.substring(i, j));
	    i=j+1;
	}
	if(i<args.length()) parm.addElement(args.substring(i, args.length()));
	String parms[] = new String[parm.size()];
	for(i=0; i<parms.length; i++) parms[i]=(String)parm.elementAt(i);
	return parms;
    }

    //----------------------------------------------------------------------------------------------------//

    final private static int OUTNUM=8;
    final private static long POWOUT=(long)Math.round(Math.pow(10, 8));
    final private static double Log10=Math.log(10);

    final private static char c[] = new char[16];
    final private static char z0='0', zm='-', zd='.', ze='e';

    /**
     * Format the double
     */

    public static String f(double d){
	char i, k;
	int log, dot;
	long aux, res;
	double abs;
	boolean end;

	if(d>0){ k=0; abs=d; }
	else if(d==0) return "0";
	else{ c[0]=zm; k=1; abs=-d; }

	log=(int)Math.floor(Math.log(abs)/Log10);
	res=(long)Math.round(abs*Math.pow(10, OUTNUM-1-log));
	if(res>=POWOUT){ res/=10; log++; }

	if(log<=-2 || log>=OUTNUM){ dot=k+1; end=true; }
	else{ dot=k+log+1; end=false; }

	k+=OUTNUM;
	for(i=0; i<=OUTNUM; i++){ if(k-i==dot) c[dot]=zd; else{ c[k-i]=(char)(z0+res%10); res/=10; } }

	while(c[k]==z0) k--;
	if(c[k]==zd) k--;
	k++;

	if(end){
	    c[k]=ze; k++;
	    if(log<0){ c[k]=zm; k++; log=-log; }
	    if(log<10) aux=1;
	    else if(log<100) aux=2;
	    else aux=3;
	    k+=aux-1;
	    for(i=0; i<aux; i++){ c[k-i]=(char)(z0+log%10); log/=10; }
	    k++;
	}

	return new String(c, 0, k);
    }

}
