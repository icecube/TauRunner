package web;
import mmc.*;
import java.io.PrintStream;
import java.applet.Applet;
import java.awt.*;
import java.awt.event.*;

/**
 * This applet class is the front-end for webmmc.
 */

public class web extends java.applet.Applet implements ActionListener{

    Propagate p;
    double vcut=1;
    double ecut=-1;
    String med="ice";
    String type="mu";

    boolean lpm=false;
    boolean cont=false;
    boolean scat=false;
    boolean time=false;
    boolean debug=true;
    boolean intr=false;

    int bs=1;
    int ph=1;
    int sh=1;
    int bb=1;
    int romb=5;
    boolean ini=false, reset=true;

    public void actionPerformed(ActionEvent event){
	double e=10, x=10, result;
	String line;

	double vcut, ecut;
	String med, type;
	boolean lpm, cont, scat, time, intr;
	int bs, ph, sh, bb;
	boolean flag, iflag;

	if(event.getSource() == panel2.button){
	    type=panel3.type.getSelectedItem();
	    if(type.equals("mpl")) type="monopole";
	    if(type.equals("all")){
		Output.DEBUG=false; panel5.debug.setState(false);
		almc lmc = new almc(panel4.text1.getText()+" -tdir="+getCodeBase()+" -user -raw "+panel4.text2.getText());
		reset=true;
		return;
	    }
	    else stdReset();

	    e=(new Double(panel4.text1.getText())).doubleValue();
	    x=(new Double(panel4.text2.getText())).doubleValue();

	    out.append("x="+Output.f(x)+" e="+Output.f(e)+" ");
	    par.append("\nplease wait ...\n");

	    lpm=panel5.lpm.getState();
	    cont=panel5.cont.getState();
	    scat=panel5.scat.getState();
	    time=panel5.time.getState();
	    debug=panel5.debug.getState();
	    intr=panel5.inter.getState();
	    switch(panel5.cuts.getSelectedIndex()){
	    case 0: vcut=1; ecut=-1; break;
	    case 1: vcut=0.05; ecut=-1; break;
	    case 2: vcut=0.01; ecut=-1; break;
	    case 3: vcut=1.e-3; ecut=-1; break;
	    default:
	    case 4: vcut=-1; ecut=500; break;
	    }
	    switch(panel5.brem.getSelectedIndex()){
	    case 0: bs=1; break;
	    case 1: bs=2; break;
	    case 2: bs=3; break;
	    default:
	    case 3: bs=4; break;
	    }
	    switch(panel5.phnu.getSelectedIndex()){
	    case 0: ph=1; sh=1; bb=3; break;
	    case 1: ph=1; sh=1; bb=4; break;
	    case 2: ph=2; sh=1; bb=3; break;
	    case 3: ph=3; sh=1; bb=1; break;
	    case 4: ph=3; sh=1; bb=2; break;
	    default:
	    case 5: ph=4; sh=2; bb=1; break;
	    }
	    med=panel3.medi.getSelectedItem();
	    switch(panel3.romb.getSelectedIndex()){
	    case 0: romb=3; break;
	    case 1: romb=4; break;
	    default:
	    case 2: romb=5; break;
	    }

	    flag=ini; ini=true; iflag=this.intr;
	    if(this.vcut!=vcut){ this.vcut=vcut; flag=false; }
	    if(this.ecut!=ecut){ this.ecut=ecut; flag=false; }
	    if(!this.med.equals(med)){ this.med=med; flag=false; }
	    if(!this.type.equals(type)){ this.type=type; flag=false; }

	    if(!flag) p = new Propagate(med, ecut, vcut, type);
	    Output.DEBUG=debug; if(!flag) iflag=false;

	    if(this.lpm!=lpm){ this.lpm=lpm; flag=false; }
	    if(this.cont!=cont){ this.cont=cont; flag=false; }
	    if(this.scat!=scat){ this.scat=scat; flag=false; }
	    if(this.time!=time){ this.time=time; flag=false; }
	    if(this.intr!=intr){ this.intr=intr; flag=false; }
	    if(this.bs!=bs){ this.bs=bs; flag=false; }
	    if(this.ph!=ph){ this.ph=ph; flag=false; }
	    if(this.sh!=sh){ this.sh=sh; flag=false; }
	    if(this.bb!=bb){ this.bb=bb; flag=false; }

	    p.sdec=true;
	    p.recc=true;
	    p.contiCorr=cont;
	    p.exactTime=time;
	    p.s.lpm=lpm;
	    p.s.b.form=bs;
	    p.s.n.form=ph;
	    p.s.n.shadow=sh;
	    p.s.n.bb=bb;
	    p.molieScat=scat;

	    Propagate.g=5;
	    if(!flag){
		if(intr){
		    err.append("\n");
		    p.interpolate("all");
		}
		else if(iflag){
		    err.append("\nSwitching off parameterizations:\n");
		    p.interpolate("none");
		}
	    }
	    Propagate.g=romb;

	    result=p.propagateTo(x*1.e2, e*1.e3);
	    line="ef="+Output.f(result>0?result*1.e-3:result*1.e-2);
	    if(time) line+=" tf="+Output.f(p.p.t);
	    out.append(line+"\n");
	    particle();
	}
    }

    void particle(){
	String outpar;
	outpar="x="+Output.f(p.p.x)+" "+"y="+Output.f(p.p.y)+" "+"z="+Output.f(p.p.z)+"\n";
	outpar+="t="+Output.f(p.p.t)+" r="+Output.f(p.p.r)+"\n";
	outpar+="theta="+Output.f(p.p.theta)+" phi="+Output.f(p.p.phi)+"\n";
	outpar+="energy="+Output.f(p.p.e)+" momentum="+Output.f(p.p.p);
	par.append(outpar+"\n");
    }

    static TextArea err;
    TextArea out, par;
    Panel1 panel1;
    Panel2 panel2;
    Panel3 panel3;
    Panel4 panel4;
    Panel5 panel5;

    public void init(){
	Output.web=true;

	String opts=getParameter("opts");
	if(opts!=null){
	    almc lmc=new almc("-tdir="+getCodeBase()+" "+opts);
	    return;
	}

	setBackground(Color.black);
	setForeground(Color.lightGray);
	setLayout(new GridLayout(3,1));

	panel1 = new Panel1(this);
	add(panel1);
	panel2 = new Panel2(this);
	add(panel2);
	err = new TextArea();
	err.setEditable(false);
	err.setBackground(Color.black);
	err.setForeground(Color.lightGray);
	add(err);
	stdReset();

	out.append(Output.version+"\n\n");
	err.append(Output.version+"\n");
	par.append(Output.version+"\n");
	Output.err.println("Running from "+getCodeBase());
    }

    private void stdReset(){
	if(reset){
	    reset=false;
	    Output.err = new PrintStream(new textStream(err));
	    Output.out = new PrintStream(new textStream(out));
	}
    }

    public String[][] getParameterInfo(){
	String pinfo[][] = {
	    {"opts", "String", "almc parameter string"},
	};
	return pinfo;
    }
}


class Panel1 extends Panel{
    Panel1(web main){
	setLayout(new GridLayout(1,2));
	main.out = new TextArea();
	main.out.setEditable(false);
	main.out.setBackground(Color.black);
	main.out.setForeground(Color.lightGray);
	add(main.out);
	main.panel3 = new Panel3(main);
	add(main.panel3);
    }
}


class Panel2 extends Panel{
    Button button;
    Panel2(web main){
	setLayout(new GridLayout(1,2));
	main.par = new TextArea();
	main.par.setEditable(false);
	main.par.setBackground(Color.black);
	main.par.setForeground(Color.lightGray);
	add(main.par);
	button = new Button("start");
	button.setBackground(Color.black);
	button.setForeground(Color.blue);
	button.setFont(new Font("Default", Font.PLAIN, 32));
	add(button);
	button.addActionListener(main);
    }
}


class Panel3 extends Panel{
    Choice type, medi, romb;
    Label Romb;
    Panel3(web main){
	main.panel4 = new Panel4(main);
	add(main.panel4);
	main.panel5 = new Panel5(main);
	add(main.panel5);
	medi = new Choice();
	medi.add("water");
	medi.add("ice");
	medi.add("salt");
	medi.add("standard rock");
	medi.add("frejus rock");
	medi.add("iron");
	medi.add("hydrogen");
	medi.add("lead");
	medi.add("uranium");
	medi.add("air");
	medi.add("mineral oil");
	medi.add("antares water");
	medi.select(1);
	add(medi);
	type = new Choice();
	type.add("mu");
	type.add("tau");
	type.add("e");
	type.add("mpl");
	type.add("all");
	add(type);
	romb = new Choice();
	romb.add("3");
	romb.add("4");
	romb.add("5");
	romb.select(2);
	add(romb);
	Romb = new Label("romb", Label.CENTER);
	add(Romb);
    }
}


class Panel4 extends Panel{
    TextField text1, text2;
    Label label1, label2;
    Panel4(web main){
	setLayout(new GridLayout(2,2));
	label1 = new Label("muon energy [GeV]", Label.CENTER);
	add(label1);
	text1 = new TextField("10000.", 5);
	text1.setBackground(Color.black);
	text1.setForeground(Color.lightGray);
	add(text1);
	label2 = new Label("propagate to [m]", Label.CENTER);
	add(label2);
	text2 = new TextField("1000.", 5);
	text2.setBackground(Color.black);
	text2.setForeground(Color.lightGray);
	add(text2);
    }
}


class Panel5 extends Panel{
    Checkbox lpm, cont, scat, time, debug, inter;
    Choice cuts, brem, phnu;
    Panel5(web main){
	setLayout(new GridLayout(3,3));
	lpm = new Checkbox("LPM");
	add(lpm);
	cont = new Checkbox("CONT");
	add(cont);
	cuts = new Choice();
	cuts.add("Vcut 1");
	cuts.add("0.05");
	cuts.add("0.01");
	cuts.add("1.e-3");
	cuts.add("E 500");
	add(cuts);
	scat = new Checkbox("SCAT");
	add(scat);
	time = new Checkbox("TIME");
	add(time);
	brem = new Choice();
	brem.add("BS: KKP");
	brem.add("BS: ABB");
	brem.add("BS: P-S");
	brem.add("BS: CSC");
	add(brem);
	debug = new Checkbox("info", true);
	add(debug);
	inter = new Checkbox("intr");
	add(inter);
	phnu = new Choice();
	phnu.add("BB 1981");
	phnu.add("BB ZEUS");
	phnu.add("BB+HARD");
	phnu.add("ALLM 91");
	phnu.add("ALLM 97");
	phnu.add("BM CKMT");
	add(phnu);
    }
}
