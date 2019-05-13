package web;
import gen.*;
import mmc.*;
import java.awt.*;
import java.awt.geom.*;
import java.io.PrintStream;
import java.util.Timer;
import java.awt.event.*;
import java.awt.AWTEvent;

/**
 * This is the graphical wrapper program for the lepton generator.
 */

public class almc extends Frame implements Runnable{

    private double R0=EarthModel.R0*1.e3+2834;
    private double R1=(EarthModel.R0+100.)*1.e3;
    private double D=1730;
    private double R=400e2;
    private double L=400e2;
    private int DNUM=2;
    private boolean f2k=false;
    private boolean all=false;

    private Graphics2D g2d;
    private AtmFlux a;
    private boolean dorun=false, ready=false;
    private boolean susp=false, texi=false;
    private int dnum=0, pnum=2;

    private Rectangle2D.Double r2d, r2p, i2d, i2p;
    private Ellipse2D.Double e2d;

    private int r2n;
    private Rectangle2D.Double[] r2e;

    private int DGN=1;
    private double xc=2;
    private long runNum, evt=0;

    private int bins=140;
    private double[] x, e, l;
    private double xmin=0.5, xmax=10;
    private double emin=0.5, emax=10;
    private double lmin=0.5, lmax=10;

    private long baseT=0;
    private Thread trun;
    private String[] args;
    private textStream myOut;

    Color cls[][]={
	{
	    new Color(70, 230, 120),
	    new Color(70, 120, 230),
	    new Color(40, 120, 70),
	    new Color(40, 70, 120),
	    new Color(100, 70, 70),
	    new Color(200, 70, 70),
	},
	{
	    new Color(150, 180, 220),
	    new Color(230, 230, 100),
	    new Color(100, 100, 255),
	    new Color(255, 100, 100),
	    new Color(230, 100, 230),
	    new Color(100, 255, 100),
	},
	{
	    new Color(140, 255, 235),
	    new Color(140, 235, 255),
	    new Color(110, 235, 140),
	    new Color(110, 140, 235),
	    new Color(210, 100, 100),
	    new Color(255, 100, 100),
	},
	{
	    new Color(40, 120, 80),
	    new Color(40, 80, 120),
	    new Color(20, 60, 40),
	    new Color(20, 40, 60),
	    new Color(50, 40, 40),
	    new Color(100, 40, 40),
	},
    };

    String str[][]={
	{
	    "mu",
	    "tau",
	    "nu_mu",
	    "nu_tau",
	    "nu_e",
	    "et al.",
	},
	{
	    "delta",
	    "brems",
	    "munu",
	    "epair",
	    "hard",
	    "et al.",
	}
    };


    public almc(String args){
	this(Output.splitString(args));
    }


    public almc(String[] args){

	for(int i=0; i<args.length; i++){
	    if(args[i].equals("-help") || args[i].equals("-h") || args[i].equals("--help")){
		Output.out.println("\n"+
"This program visualizes atmospheric lepton fluxes at 0.1-4000 TeV\n"+
"                       -f2k  output f2k lines of generated events\n"+
"                       -all  display all generated events\n"+
"                       -dnum=[num. of events kept before refresh]\n"+
"                       -dgn=[N: 2N+1 sets the max number of  and]\n"+
"                       -xc=[sets distance between  spread tracks]\n"+
"atmflux help lines follow:");
		new AtmFlux().initAtmFlux(args, true);
		return;
	    }
	    else if(args[i].equals("-f2k")){
		f2k=true;
	    }
	    else if(args[i].equals("-all")){
		all=true;
	    }
	    else if(args[i].startsWith("-dnum=")){
		try{
		    DNUM=(int)Double.parseDouble(args[i].substring(6));
		}catch(Exception error){
		    DNUM=2;
		}
	    }
	    else if(args[i].startsWith("-dgn=")){
		try{
		    DGN=(int)Double.parseDouble(args[i].substring(5));
		}catch(Exception error){
		    DGN=1;
		}
	    }
	    else if(args[i].startsWith("-xc=")){
		try{
		    xc=Double.parseDouble(args[i].substring(4));
		}catch(Exception error){
		    xc=2;
		}
	    }
	}

	if(DNUM<=0) DNUM=2;
	if(DGN<0 || DGN>10) DGN=1;
	if(xc<0 || xc>10) xc=2;

	setSize(600, 600);
	setResizable(false);

	String title="all lepton propagation with mmc";
	setTitle(title);

	if(Output.web) Output.texi=false;
	else{
	    Panel head = new Panel();
	    head.add(new Label(title));
	    head.setForeground(Color.red);
	    add(head, BorderLayout.NORTH);
	}

	this.args=args;

	x = new double[bins];
	for(int i=0; i<bins; i++) x[i]=0;
	e = new double[bins];
	for(int i=0; i<bins; i++) e[i]=0;
	l = new double[bins];
	for(int i=0; i<bins; i++) l[i]=0;

	e2d = new Ellipse2D.Double(50, 60, 500, 500);
	r2d = new Rectangle2D.Double(400, 60, 160, 160);
	r2p = new Rectangle2D.Double(350, 410, 220, 150);
	i2d = new Rectangle2D.Double(401, 61, 159, 159);
	i2p = new Rectangle2D.Double(351, 411, 219, 149);

	r2n=3;
	r2e = new Rectangle2D.Double[r2n];
	r2e[0] = new Rectangle2D.Double(0, 410, 350, 190);
	r2e[1] = new Rectangle2D.Double(0, 221, 600, 189);
	r2e[2] = new Rectangle2D.Double(0, 0, 400, 221);

	Output.err = new PrintStream(myOut = new textStream(this));
	if(!f2k) Output.out=Output.err;

	setVisible(true);

	enableEvents(AWTEvent.MOUSE_EVENT_MASK | AWTEvent.WINDOW_EVENT_MASK);

	trun = new Thread(this);
	trun.start();
    }


    protected void processMouseEvent(MouseEvent e){
	if(e.MOUSE_CLICKED==e.getID()){
	    double x=e.getX(), y=e.getY();
	    if(x>=400 && x<=560 && y>=60 && y<=220){
		susp=!susp;
	    }
	    else if(x>=350 && x<=570 && y>=410 && y<=560){
		if(ready) runpaint(3, null);
	    }
	    else if((x-300)*(x-300)+(y-310)*(y-310)<250*250){
		all=!all;
		if(!all) repaint();
	    }
	}
	super.processMouseEvent(e);
    }


    protected void processWindowEvent(WindowEvent e){
	if(e.WINDOW_CLOSING==e.getID()){
	    texi=true;
	    Output.texi=true;
	    disableEvents(AWTEvent.MOUSE_EVENT_MASK | AWTEvent.WINDOW_EVENT_MASK);
	    dispose();
	    try{
		System.exit(0);
	    } catch(Exception ex) {
	    }
	}
	super.processWindowEvent(e);
    }


    public static void main(String[] args){
	new almc(args);
    }


    private void xplot(int i, int q){
	double f;
	switch(q){
	case 0:
	    f=x[i];
	    f=f>0?Math.log(f/xmin)/Math.log(xmax/xmin):0;
	    g2d.setColor(Color.green);
	    break;
	case 1:
	    f=e[i];
	    f=f>0?Math.log(f/emin)/Math.log(emax/emin):0;
	    g2d.setColor(Color.orange);
	    break;
	case 2:
	    f=l[i];
	    f=f>0?Math.log(f/emin)/Math.log(lmax/lmin):0;
	    g2d.setColor(Color.magenta);
	    break;
	default: f=0;
	}
	g2d.setStroke(new BasicStroke(1f));
	g2d.draw(new Line2D.Double(px((i+0.)/bins), py(f), px((i+1.)/bins), py(f)));
    }


    private void xreplot(int q){
	g2d.setClip(i2p);
	g2d.setColor(Color.black);
	g2d.fill(i2p);

	for(int i=0; i<bins; i++) xplot(i, q);
	g2d.draw(new Line2D.Double(px(0), py(0), px(1), py(0)));
	switch(q){
	case 0:
	    g2d.drawString("zenith angle 0 - 180 degrees", (int)px(0.03), (int)py(-0.1));
	    for(int n=0; n<=18; n++) g2d.draw(new Line2D.Double(px(n/18.), py(0.01), px(n/18.), py(-0.02)));
	    break;
	case 1:
	    g2d.drawString("energy deposited .1 GeV - 1 PeV", (int)px(0.01), (int)py(-0.1));
	    for(int n=0; n<=7; n++) g2d.draw(new Line2D.Double(px(n/7.), py(0.01), px(n/7.), py(-0.02)));
	    break;
	case 2:
	    g2d.drawString("initial energy 1 GeV - 1 EeV", (int)px(0.05), (int)py(-0.1));
	    for(int n=0; n<=9; n++) g2d.draw(new Line2D.Double(px(n/9.), py(0.01), px(n/9.), py(-0.02)));
	    break;
	}
    }


    private void aplot(int i, int q){
	g2d.setClip(i2p);
	switch(q){
	case 0:
	    if(x[i]>=xmax){
		xmax=10*xmax;
		xreplot(q);
	    }
	    else xplot(i, q);
	    break;
	case 1:
	    if(e[i]>=emax){
		emax=10*emax;
		xreplot(q);
	    }
	    else xplot(i, q);
	    break;
	case 2:
	    if(l[i]>=lmax){
		lmax=10*lmax;
		xreplot(q);
	    }
	    else xplot(i, q);
	    break;
	}
    }


    private boolean runchecks(){
	if(force){
	    g2d = (Graphics2D)getGraphics();
	    incnum(false);
	    force=false;
	}
	if(texi) if(f2k){
	    Output.out.println("TEND ? ? ?");
	    Output.out.println("END");
	}
	return texi;
    }


    public void run(){
	a = new AtmFlux();
	a.setup(args);
	D=a.D;
	R=a.radius*1.e2;
	L=a.length*1.e2/2;
	if(texi) return;
	ready=true;
	repaint();
	if(f2k){
	    Output.out.println("V 2000.1.2");
	    Output.out.println("HI  almc: event display tool");
	    Output.out.println("USER_DEF almc_e Name Eini Elost Theta");
	    Output.out.println("TBEGIN ? ? ?");
	}

	while(true){
	    runNum=0;
	    while(runNum<100000){
		while(susp) try{
		    if(runchecks()) return;
		    trun.sleep(100);
		} catch(Exception ex) {
		}
		if(runchecks()) return;
		runpaint(1, null);
	    }
	    repaint();
	}
    }


    public void paint(Graphics gr){
	runpaint(2, gr);
    }


    private synchronized void runpaint(int i, Graphics gr){
	switch(i){
	case 1: arun(); break;
	case 2: mypaint(gr); break;
	case 3: incnum(true); break;
	}
    }


    private boolean force=false;

    private void incnum(boolean force){
	long t=System.currentTimeMillis();
	if(force){
	    baseT=t-10001;
	    this.force=true;
	}
	else if(t-baseT>10000){
	    baseT=t;
	    if(pnum>=2) pnum=0; else pnum++;
	    xreplot(pnum);
	}
    }


    private void arun(){
	dorun=false;
	g2d = (Graphics2D)getGraphics();

	incnum(false);

	Particle[] pSet=a.createNext();
	if(a.elost>0){
	    int i=0, k;

	    k=Math.min((int)Math.floor(bins*(180-pSet[0].theta)/180), bins-1);
	    x[k]++;
	    if(pnum == 0) i=k;

	    k=Math.min((int)Math.floor(bins*Math.min(Math.max(Math.log(a.elost/.1)/Math.log(1.e7), 0), 1)), bins-1);
	    e[k]++;
	    if(pnum == 1) i=k;

	    k=Math.min((int)Math.floor(bins*Math.min(Math.max(Math.log(pSet[0].e/1.e3)/Math.log(1.e9), 0), 1)), bins-1);
	    l[k]++;
	    if(pnum == 2) i=k;

	    aplot(i, pnum);

	    if(dnum>=DNUM){
		dnum=0;
		g2d.setClip(i2d);
		g2d.setColor(Color.black);
		g2d.fill(i2d);
	    }
	    dnum++;
	}
	if((all && a.elost>=0) || (!all && a.elost>0)){
	    runNum++;
	    if(f2k) Output.out.println("EM "+evt+" 1 1970 0 "+runNum+" 0");
	    g2d.setStroke(new BasicStroke(1f));

	    int dgr=-1, dgs=-1;
	    for(int i=0; i<pSet.length; i++){
		Particle p=pSet[i];
		int type;
		Line2D line;
		double gx, gy;
		String name=p.name;

		if(name.equals("mu-") || name.equals("mu+")) type=0;
		else if(name.equals("tau-") || name.equals("tau+")) type=1;
		else if(name.equals("nu_mu") || name.equals("~nu_mu")) type=2;
		else if(name.equals("nu_tau") || name.equals("~nu_tau")) type=3;
		else if(name.equals("nu_e") || name.equals("~nu_e")) type=4;
		else type=5;

		if(p.r>0){
		    double dr=1.e7;
		    double xi, xj, xf, yi, yj, yf, aux1, aux2;

		    gx=xc*p.costh; gy=xc*p.sinth;

		    xi=Math.sqrt(p.x*p.x+p.y*p.y);
		    if(p.x<0) xi=-xi;
		    yi=p.z;
		    aux1=p.x+p.r*p.sinth*p.cosph;
		    aux2=p.y+p.r*p.sinth*p.sinph;
		    xf=Math.sqrt(aux1*aux1+aux2*aux2);
		    if(aux1<0) xf=-xf;
		    yf=yi+p.r*p.costh;

		    if((xi-xf)*(yi-yf)*gx*gy<0) gx*=-1;

		    if(p.r>2*dr){
			if(dgr>=DGN) dgr=-DGN; else dgr++;
		    }
		    else{
			gx=0; gy=0;
		    }

		    if(p.r<dr) dr=p.r/2;
		    aux1=p.x+(p.r-dr)*p.sinth*p.cosph;
		    aux2=p.y+(p.r-dr)*p.sinth*p.sinph;
		    xj=Math.sqrt(aux1*aux1+aux2*aux2);
		    if(aux1<0) xj=-xj;
		    yj=yi+(p.r-dr)*p.costh;

		    for(int j=0; j<r2n; j++){
			g2d.setClip(r2e[j]);

			line = new Line2D.Double(gx(xi), gy(yi), gx(xf)+dgr*gx, gy(yf)+dgr*gy);
			g2d.setColor(cls[0][type]);
			g2d.draw(line);

			line = new Line2D.Double(gx(xj)+dgr*gx*(1-dr/p.r), gy(yj)+dgr*gy*(1-dr/p.r), gx(xf)+dgr*gx, gy(yf)+dgr*gy);
			g2d.setColor(cls[2][type]);
			g2d.draw(line);
		    }
		}

		if(a.elost>0){
		    g2d.setClip(i2d);
		    if(p.r==0){
			if(name.equals("delta")) type=0;
			else if(name.equals("brems")) type=1;
			else if(name.equals("munu")) type=2;
			else if(name.equals("epair")) type=3;
			else if(name.equals("hadr")) type=4;
			else type=5;
			g2d.setColor(cls[1][type]);
			double r=Math.min(Math.max(5*Math.log(p.e/500)/Math.log(1.e5/500)+2, 2), 7);
			Ellipse2D e2d = new Ellipse2D.Double(dx(p.x)-r/2, dy(p.z)-r/2, r, r);
			g2d.fill(e2d);
		    }
		    else{
			boolean flag=true;
			double x1, x2, y1, y2, x1f, x2f, y1f, y2f, aux;

			x1=p.x;
			x2=p.x+p.r*p.sinth*p.cosph;
			y1=p.z;
			y2=p.z+p.r*p.costh;

			// this horror is necessary since lines fail to plot correctly with large coordinates
			if(Math.max(Math.abs(x2),Math.abs(x1))>Math.max(Math.abs(y2),Math.abs(y1))){

			    if(x1<-R){
				if(x2>-R){
				    aux=(-R-x1)/(x2-x1);
				    x1+=aux*(x2-x1);
				    y1+=aux*(y2-y1);
				    if(x2>R){
					aux=(R-x2)/(x1-x2);
					x2+=aux*(x1-x2);
					y2+=aux*(y1-y2);
				    }
				}
				else flag=false;
			    }
			    else if(x1>R){
				if(x2<R){
				    aux=(R-x1)/(x2-x1);
				    x1+=aux*(x2-x1);
				    y1+=aux*(y2-y1);
				    if(x2<-R){
					aux=(-R-x2)/(x1-x2);
					x2+=aux*(x1-x2);
					y2+=aux*(y1-y2);
				    }
				}
				else flag=false;
			    }
			    else{
				if(x2>R){
				    aux=(R-x2)/(x1-x2);
				    x2+=aux*(x1-x2);
				    y2+=aux*(y1-y2);
				}
				else if(x2<-R){
				    aux=(-R-x2)/(x1-x2);
				    x2+=aux*(x1-x2);
				    y2+=aux*(y1-y2);
				}
			    }

			}
			else{

			    if(y1<-L){
				if(y2>-L){
				    aux=(-L-y1)/(y2-y1);
				    x1+=aux*(x2-x1);
				    y1+=aux*(y2-y1);
				    if(y2>L){
					aux=(L-y2)/(y1-y2);
					x2+=aux*(x1-x2);
					y2+=aux*(y1-y2);
				    }
				}
				else flag=false;
			    }
			    else if(y1>L){
				if(y2<L){
				    aux=(L-y1)/(y2-y1);
				    x1+=aux*(x2-x1);
				    y1+=aux*(y2-y1);
				    if(y2<-L){
					aux=(-L-y2)/(y1-y2);
					x2+=aux*(x1-x2);
					y2+=aux*(y1-y2);
				    }
				}
				else flag=false;
			    }
			    else{
				if(y2>R){
				    aux=(R-y2)/(y1-y2);
				    x2+=aux*(x1-x2);
				    y2+=aux*(y1-y2);
				}
				else if(y2<-R){
				    aux=(-R-y2)/(y1-y2);
				    x2+=aux*(x1-x2);
				    y2+=aux*(y1-y2);
				}
			    }

			}

			if(flag){
			    if(dgs>=DGN) dgs=-DGN; else dgs++;
			    aux=p.sinth*p.sinph;
			    aux=Math.sqrt(1-aux*aux);
			    if(aux>0){
				gx=xc*p.costh/aux; gy=xc*p.sinth*p.cosph/aux;
			    }
			    else{
				gx=0; gy=0;
			    }
			    line = new Line2D.Double(dx(x1), dy(y1), dx(x2)+dgs*gx, dy(y2)+dgs*gy);
			    g2d.setColor(cls[3][type]);
			    g2d.draw(line);
			}
		    }
		}

		if(f2k) Output.out.println("TR "+p.gens+" "+p.igen+" "+p.name+" "+
					   Output.f(p.x*1.e-2)+" "+Output.f(p.y*1.e-2)+" "+Output.f(p.z*1.e-2)+" "+
					   Output.f(180-p.theta)+" "+Output.f(p.phi<180?p.phi+180:p.phi-180)+" "+
					   Output.f(p.r*1.e-2)+" "+Output.f(p.e*1.e-3)+" "+Output.f(p.t*1.e9));
	    }
	    if(f2k) Output.out.println("US almc_e "+a.Name+" "+Output.f(a.eini)+" "+
				       Output.f(a.elost)+" "+Output.f(a.thini)+"\nEE");
	}
	dorun=true;
    }


    private double gx(double x){
	return 300+250*x/(R1*1.e2);
    }


    private double gy(double y){
	return 310-250*(y+(R0-D)*1.e2)/(R1*1.e2);
    }


    private double dx(double x){
	return 480+80*x/(R);
    }


    private double dy(double y){
	return 140-80*y/(L);
    }


    private double px(double x){
	return 355+210*x;
    }


    private double py(double y){
	return 545-130*y;
    }


    private void mypaint(Graphics gr){
	g2d = (Graphics2D)gr;

	setBackground(Color.black);
	setForeground(Color.lightGray);

	g2d.setClip(null);

	myOut.reflush();
	if(ready){
	    g2d.setStroke(new BasicStroke(2f));
	    for(int i=0; i<250; i++){
		double rho=a.eM.rho(i*(R1*1.e-3)/250);
		int r, g, b;
		if(rho>1){
		    r=Math.min(Math.max((int)Math.floor(120-70*Math.log(rho/1.)/Math.log(14/1.)), 0), 255);
		    g=r/2;
		    b=r/2;
		}
		else{
		    r=0;
		    g=0;
		    b=255;
		}
		g2d.setColor(new Color(r, g, b));
		g2d.draw(new Ellipse2D.Double(300-i, 310-i, 2*i, 2*i));
	    }
	}
	else{
	    g2d.setStroke(new BasicStroke(4f));
	    g2d.setColor(new Color(0, 0, 230));
	    g2d.draw(e2d);
	}

	g2d.setColor(Color.black);
	g2d.fill(r2d);
	g2d.setStroke(new BasicStroke(1f));
	g2d.setColor(Color.blue);
	g2d.draw(r2d);

	g2d.setColor(Color.black);
	g2d.fill(r2p);
	g2d.setStroke(new BasicStroke(1f));
	g2d.setColor(Color.blue);
	g2d.draw(r2p);

	for(int i=0; i<6; i++){
	    g2d.setColor(cls[0][i]);
	    g2d.drawString(str[0][i], 30, 460+i*20);
	    g2d.setColor(cls[1][i]);
	    g2d.drawString(str[1][i], 30, 70+i*20);
	}

	if(ready) xreplot(pnum);

	else g2d.drawString("initializing, please wait ...", (int)px(0.05), (int)py(-0.1));
	runNum=0;
    }

}
